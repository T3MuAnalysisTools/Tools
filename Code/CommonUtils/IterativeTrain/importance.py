import pandas as pd
import matplotlib.pyplot as plt

# List of file names
file_names = ['tauha', 'tauhb', 'taumu', 'taue']
title_type = ["{h,A}$","{h,B}$","{\mu}$","{e}$"]

# Plotting
for i, file_name in enumerate(file_names):
    # Load data from data.dat
    data = pd.read_csv('importances/'+file_name+'.txt', engine='python', header=0, delim_whitespace=True)

    # Plot horizontal bar chart
    plt.figure(figsize=(10, 8))
    plt.barh(data['Variable'], data['Value'], color='skyblue')
    plt.xlabel('Variable Importance')
    plt.title(r'Category $\tau_'+title_type[i])
    plt.gca().invert_yaxis()  # Invert y-axis to have highest importance at the top
    
    # Adjust margins to prevent cutting off variable names
    plt.subplots_adjust(left=0.25)
    
    plt.savefig('importances/'+file_name+'.png', dpi=300)  # Save each plot with a different name
    plt.show()
