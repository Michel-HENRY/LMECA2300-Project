import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

filename = "results.csv"
data = np.loadtxt(filename, delimiter=',');
print(np.shape(data));

iter = data[0,0]
x = data[:,1:3]
u = data[:,3:5]
uBoundary = 5

r = la.norm(x,axis=1)
ur_num = la.norm(u,axis = 1)
ur_the = uBoundary*r

fig, ax = plt.subplots()
ax.plot(r,ur_num,'b--', lineWidth = 3)
ax.plot(r,ur_the,'r:' ,lineWidth = 3)
ax.set_xlabel('Radius')
ax.set_ylabel('Radial Velocity')
ax.set_xlim([0,1])
ax.set_ylim([0,5])
ax.legend(["Numerical Value","Theoretical velocity"])
plt.title('Radial Velocity in a cylinder')

plt.show()
