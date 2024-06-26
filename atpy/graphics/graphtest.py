import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2)
axes[0][0].plot([1, 2, 3], [4, 5, 6])
axes[0][0].axis('off')
axes[0][1].plot([1, 2, 3], [4, 5, 6])
axes[0][1].axis('off')
axes[1][0].plot([1, 2, 3], [4, 5, 6])
axes[1][0].axis('off')
axes[1][1].plot([1, 2, 3], [4, 5, 6])
# axes[1][1].xticks().set_visible(False)
axes[1][1].spines["bottom"].set_visible(False)
axes[1][1].set_xlabel("x[m]")
plt.show()
