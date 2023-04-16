#%%
import subprocess
#%%
# original file
data_folder="/home/vipink/Documents/FANTOM6/data/RADICL_data/Neuron/replicate1/raw/"
radicl_file = "15.GGCTAC.2.bed"
#%%
f_out = open(f"{data_folder}RADICL_intra.bed","w")
subprocess.run(['awk',
                f"""$1 == $7 """,
                f'{data_folder}{radicl_file}'],shell=False,stdout=f_out)
# %%
