"""
The DeepCaller network.

One recurrent encoder is shared by the `ploidy` allele slots of a locus. The
slot representations then attend to each other, and an autoregressive decoder
walks the slots in order, emitting a copy-number distribution for each one while
keeping track of how many copies are still unassigned.

TensorFlow is imported at module level, so this module must only be imported
after any process pool has been created: forking a process that already carries
TensorFlow's thread pools is not safe.
"""

import numpy as np
import tensorflow as tf
from tensorflow.keras import initializers, layers, regularizers

from .config import FEATURE_DIM


class DeepCallerModel(tf.keras.Model):
    """Copy-number predictor over the per-allele pileup encodings of one locus."""

    def __init__(self, ploidy, seed=42, l2_strength=0.01, num_attn_heads=1,
                 lstm_cells=64, **kwargs):
        super().__init__(**kwargs)
        self.ploidy = ploidy
        self.lstm_out_dim = lstm_cells

        # Keras defaults (glorot_uniform + orthogonal) are kept for the LSTM;
        # He initialisation is only appropriate for the ReLU layers below.
        self.LSTM = layers.LSTM(lstm_cells, return_sequences=False)

        # Cross-slot self-attention: every allele slot gets to see the others
        # before its copy number is decided.
        self.cross_head_attn = layers.MultiHeadAttention(
            num_heads=num_attn_heads, key_dim=self.lstm_out_dim // num_attn_heads)
        self.attn_norm = layers.LayerNormalization()

        # Autoregressive decoder, shared across the alt1 .. altP slots.
        self.remaining_embed = layers.Dense(
            16, activation="relu",
            kernel_initializer=initializers.HeNormal(seed=seed),
            kernel_regularizer=regularizers.l2(l2_strength))
        self.Dense_copy = layers.Dense(
            32, activation="relu",
            kernel_initializer=initializers.HeNormal(seed=seed),
            kernel_regularizer=regularizers.l2(l2_strength))
        self.Output_copy = layers.Dense(
            self.ploidy + 1,
            kernel_initializer=initializers.HeNormal(seed=seed),
            dtype="float32")

    def call(self, x, training=False):
        """
        Args:
            x: (batch, ploidy, window, features) float32 tensor.

        Returns:
            probs: (batch, ploidy, ploidy + 1) copy-number distributions.
            head_mask: (batch, ploidy) bool, False where a slot carries no reads.
        """
        batch_size = tf.shape(x)[0]
        n_groups = self.ploidy

        # Every slot is encoded independently by the shared LSTM.
        x_flat = tf.reshape(x, (-1, tf.shape(x)[2], tf.shape(x)[3]))
        h = self.LSTM(x_flat)
        h = tf.reshape(h, (batch_size, n_groups, self.lstm_out_dim))

        # An all-zero slot means "no read supported this allele"; it is masked
        # out of the attention and never consumes a copy.
        head_mask = tf.reduce_any(tf.not_equal(x, 0), axis=[2, 3])
        h = tf.where(head_mask[:, :, tf.newaxis], h, tf.zeros_like(h))

        attn_mask = head_mask[:, tf.newaxis, :]
        h_attn = self.cross_head_attn(
            query=h, key=h, value=h, attention_mask=attn_mask, training=False)
        h_ctx = self.attn_norm(h + h_attn)

        # Decode alt1 .. altP, carrying the number of unassigned copies forward.
        remaining = tf.fill((batch_size,), tf.cast(self.ploidy, tf.float32))
        probs_list = []
        copy_values = tf.range(self.ploidy + 1, dtype=tf.float32)

        for step in range(self.ploidy):
            head_repr = h_ctx[:, step, :]
            remaining_feat = self.remaining_embed(remaining[:, tf.newaxis])
            z = self.Dense_copy(tf.concat([head_repr, remaining_feat], axis=-1))
            logits = self.Output_copy(z)

            # Copy numbers above the remaining budget are impossible.
            valid = tf.cast(copy_values[tf.newaxis, :] <= remaining[:, tf.newaxis],
                            tf.float32)
            logits = tf.where(valid > 0, logits, tf.fill(tf.shape(logits), -1e9))
            probs = tf.nn.softmax(logits, axis=-1)
            probs_list.append(probs)

            # Inference: the budget is updated from the model's own argmax,
            # there is no teacher forcing here.
            is_masked = tf.logical_not(head_mask[:, step])
            step_dosage = tf.cast(tf.argmax(probs, axis=-1), tf.float32)
            step_dosage = tf.where(is_masked, tf.zeros_like(step_dosage), step_dosage)
            remaining = remaining - step_dosage

        return tf.stack(probs_list, axis=1), head_mask

    def predict_in_batches(self, features, batch_size):
        """
        Run inference over a large encoding array in fixed-size batches.

        Each batch is cast to float32 on the fly, so the caller can keep the
        encodings in float16 and halve their resident size.
        """
        n = features.shape[0]
        if n == 0:
            return (np.empty((0, self.ploidy, self.ploidy + 1), dtype=np.float32),
                    np.empty((0, self.ploidy), dtype=bool))

        probs_chunks, mask_chunks = [], []
        for start in range(0, n, batch_size):
            batch = tf.convert_to_tensor(
                features[start:start + batch_size], dtype=tf.float32)
            probs, mask = self(batch, training=False)
            probs_chunks.append(probs.numpy())
            mask_chunks.append(mask.numpy())

        return (np.concatenate(probs_chunks, axis=0),
                np.concatenate(mask_chunks, axis=0))


def load_model(weights_path, window_length, ploidy):
    """
    Build the graph and restore the weights.

    A dummy forward pass is required first: the variables are only created once
    the layers have seen an input shape.
    """
    model = DeepCallerModel(ploidy=ploidy)
    model(tf.zeros((1, ploidy, window_length, FEATURE_DIM), dtype=tf.float32),
          training=False)
    model.load_weights(weights_path)
    return model


def configure_threading(num_threads):
    """Pin TensorFlow to the CPU and bound its intra/inter op parallelism."""
    tf.config.set_visible_devices([], "GPU")
    tf.config.threading.set_intra_op_parallelism_threads(num_threads)
    tf.config.threading.set_inter_op_parallelism_threads(num_threads)
