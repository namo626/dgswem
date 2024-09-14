import numpy as np
import pandas as pd
import matplotlib.tri as tri
import os


def read_triangulation(fname: str) -> (tri.Triangulation, np.ndarray):
    """ Read a fort.14 file and return its triangulation data as a
    matplotlib Triangulation object.
    """
    df = pd.read_csv(fname, sep="\\s+", names=list('abcde'),on_bad_lines='skip')

    # Extract the triangulation, and drop them from the main dataframe
    triangles = df.loc[df['e'].notnull()].iloc[:,2:].values
    # Minus one due to Python's indexing
    triangles = triangles.astype(int) - 1

    df = df[~df['e'].notnull()]

    # Extract the nodal coordinates
    nodes = df.loc[df['d'].notnull()][['b', 'c']].values.astype(float)

    # Extract the bathymetry
    bathy = df.loc[df['d'].notnull()]['d'].values.astype(float)


    tr = tri.Triangulation(nodes[:,0], nodes[:,1], triangles)
    return tr, bathy


def read_63(fname: str, nodes: list, total_nodes: int, col: int=1) -> list:
    """ Given a fort.63-like file, extract the nodal values of the given list of nodes
    and return them as a list. Values from different timestamps are concatenated without
    space, i.e. the last node at timestep n is immediately followed by first node
    at timestep n+1. If the file is a fort.64 and we want the y-velocity values,
    use col=2.
    """
    df = pd.read_csv(fname, skiprows=2, sep="\\s+", index_col=False, header=None, names=list('abc'))

    # drop the timestamp rows
    ts = list(filter(lambda x: (x % (total_nodes+1) == 0), range(df.shape[0])))
    df = df.drop(ts)

    df.iloc[:,0] = df.iloc[:,0].astype(int)

    # Filter the nodes as specified
    cross_elev = df.loc[df.iloc[:,0].isin(nodes)]
    cross_elev.iloc[:,0] = cross_elev.iloc[:,0].astype(int)

    return cross_elev.iloc[:,col].values

def read_63_all(fname: str, meshname: str):
    tr,_ = read_triangulation(meshname)
    total_nodes = len(tr.x)
    node_list = list(range(1, total_nodes+1))

    d1 = read_63(fname, nodes=node_list, total_nodes=total_nodes)
    return tr, d1
