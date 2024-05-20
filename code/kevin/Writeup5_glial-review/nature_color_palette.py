import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Define the color palette
donor_col_palette = {
    "D:10": mcolors.to_hex((95/255, 227/255, 255/255)),
    "D:13": mcolors.to_hex((3/255, 170/255, 255/255)),
    "D:16": mcolors.to_hex((93/255, 157/255, 255/255)),
    "D:1": mcolors.to_hex((0/255, 62/255, 161/255)),
    "D:20": mcolors.to_hex((5/255, 75/255, 255/255)),
    "D:17": mcolors.to_hex((109/255, 157/255, 255/255)),
    "D:22": mcolors.to_hex((0/255, 22/255, 170/255)),
    # break, cases below
    "D:11": mcolors.to_hex((255/255, 102/255, 162/255)),
    "D:2": mcolors.to_hex((255/255, 56/255, 156/255)),
    "D:3": mcolors.to_hex((191/255, 0/255, 86/255)),
    "D:12": mcolors.to_hex((255/255, 46/255, 118/255)),
    "D:4": mcolors.to_hex((255/255, 102/255, 178/255)),
    "D:5": mcolors.to_hex((255/255, 86/255, 137/255)),
    "D:8": mcolors.to_hex((174/255, 0/255, 60/255))
}

batch_col_palette = {
    "1": mcolors.to_hex((255/255, 135/255, 73/255)),
    "2": mcolors.to_hex((205/255, 175/255, 255/255))
}