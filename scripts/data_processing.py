import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('../build/data.csv')
# print(df)

fig, axes = plt.subplots(nrows=2, ncols=2)

df.plot(ax = axes[0, 0], grid = True, x = 't', y = ['ex', 'ey', 'ez'])

df.plot(ax = axes[0, 1], grid = True, x = 't', y = ['ephi', 'etheta', 'epsi'])

df.plot(ax = axes[1, 0], grid = True, x = 't', y = ['fz'])

df.plot(ax = axes[1, 1], grid = True, x = 't', y = ['z'])

plt.show()