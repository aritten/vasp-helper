#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np

data = np.array([11.15,7.29])

fig, ax = plt.subplots(1, 1, figsize=(6,6))
ax.bar(['vasphelper', 'Old Script'], data, width = .4, color=['blue','green'], edgecolor='black',linewidth=2)
ax.set_ylabel('Time (sec)')
ax.set_xlabel('Method')

# sub_ax = ax.inset_axes(
#     bounds=(0.15,0.4,0.35,0.5)
# )

# # sub_ax.bar(['vasphelper'], data[0], width = 0.3, color='skyblue', edgecolor='black', linewidth=2)
# # sub_ax.set_ylabel('Time (sec)')
# # sub_ax.set_xlim(-0.3, 0.3)

plt.show()