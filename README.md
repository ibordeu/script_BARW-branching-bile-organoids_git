# Human branching cholangiocyte organoids recapitulate functional bile duct formation. 
## Roos et al. Cell Stem Cell (2022)

Simulation of the three-dimensional (3d) branching-annihilating random walk
considering: tip-duct repulsion, a refractory period for branching and
terminations, and a typical tip size at branching. For details, see
Methods in paper.
 
### INPUTS: see MODEL PARAMETERS below, and Methods. The default parameters 
provided, correspond to the parameters used to simulate the BRCO condition. 

### OUTPUTS:
edge_list: Array of size number_of_links*3, containing the source node id (column 1), target node id (column 2), and distance from source to target (column 3). For all nodes in the simulation.

node_positions : Array of size number_of_nodes*4. Column 1: node id. Columns 2-4: coordinates (x,y,z) for each node, respectively.  

edge_list_tree: Reduced edge_list, where each node corresponds to either a branching point or a termination point or the root node

node_positions_tree : node positions for the reduced branching tree.
