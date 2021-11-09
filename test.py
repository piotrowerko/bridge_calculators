import numpy as np


# 4x + 3y + 2z = 25
# -2x + 2y + 3z = -10
# 3x -5y + 2z = -4

A = np.array([[4, 3, 2], [-2, 2, 3], [3, -5, 2]])
B = np.array([25, -10, -4])
X = np.linalg.inv(A).dot(B)

print(X)