# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
# %%
data=pd.read_csv("Trial_file.csv")

# %%
sns.displot(data['Seq'], kde=True)

# %%
