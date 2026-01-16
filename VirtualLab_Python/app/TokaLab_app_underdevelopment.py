# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 17:53:36 2026

@author: ricca
"""

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

st.title("La mia prima app")

x = st.slider("Scegli un valore", 0.0, 10.0, 5.0)

t = np.linspace(0, 10, 200)
y = np.sin(x * t)

fig, ax = plt.subplots()
ax.plot(t, y)
st.pyplot(fig)