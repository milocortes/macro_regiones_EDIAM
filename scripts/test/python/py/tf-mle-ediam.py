import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow_probability import distributions as tfd
from scipy.integrate import odeint

NUM_SAMPLES = 1000
device = 'cpu:0' # These experiments do not require the GPU. Normally, 'gpu:0' if tf.test.is_gpu_available() else 'cpu:0' should be used.
initial_states = tf.convert_to_tensor([1., 1.], dtype=tf.float64); initial_states
t = tf.linspace(0., 1., num=NUM_SAMPLES); t.shape


class Lambda(tf.keras.Model):

  def call(self, t, y):
    # y now reprents the a vector of [u, v]
    u, v = y[0], y[1]

    du_dt = v
    dv_dt = 5 * v - 6 * u

    return tf.stack([du_dt, dv_dt])  # vector of shape [2]
# https://rlhick.people.wm.edu/posts/custom-likes-tensorflow.html
# https://mauriciotejada.com/programacionjulia/aplicaciones-iii.html#estimaci%C3%B3n-por-m%C3%A1ximo-verosimilitud
# https://jeffpollock9.github.io/maximum-likelihood-estimation-with-tensorflow-probability-and-stan-take-2/
