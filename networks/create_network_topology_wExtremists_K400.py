import numpy as np
import logging
import random as rg
import sys
from timeit import default_timer as timer
import csv
import os


def is_even(x):
    # Input must be an int
    # Output is a boolean
    if type(x) is not int:
        return None

    y = np.dtype(bool)

    if x % 2 == 0:
        y = True
    else:
        y = False

    return y


def export_world_as_CSV(W, Nm, Nx1, Nx2, file_name):
    with open(file_name, mode='w', newline='') as file:
        Wn = W.shape[0]
        header = ['']
        all_rows = []
        file_writer = csv.writer(file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        # Construct the header
        for v in range(0, Wn):
            if v < Nm:
                header.append(str("m") + str(v))
            elif v < Nm + Nx1:
                header.append(str("px") + str(v - Nm))
            else:
                header.append(str("nx") + str(v - Nm - Nx1))

        # Construct a list of rows
        for v in range(0, Wn):
            if v < Nm:
                row = [str("m") + str(v)]
            elif v < Nm + Nx1:
                row = [str("px") + str(v - Nm)]
            else:
                row = [str("nx") + str(v - Nm - Nx1)]
            for n in range(0, Wn):
                row.append(W[v][n])
            all_rows.append(row)

        # Write to file
        file_writer.writerow(header)
        file_writer.writerows(all_rows)
    return


def get_ring_lattice(Nv, Kv):
    """
    Create a ring lattice of Nv nodes, with degree Kv.

    :param Nv: Number of vertices. int in [0,1000].
    :param Kv: Degree of a vertex. int in [0,Nv-1]; must be an even number.
    :return: N: Matrix as the ring-lattice network
    """
    if type(Nv) is not int:
        logging.error("Error: Nv should be an int data type.")
        return None

    if type(Kv) is not int:
        logging.error("Error: Kv should be an int data type.")
        return None

    if Nv < 0 or Nv > 1000:
        logging.error("Error: Nv should be in [0,1000].")
        return None

    if Kv < 0 or Kv > Nv - 1:
        logging.error("Error: Kv should be in [0, Nv-1].")
        return None

    if is_even(Kv) is not True:
        logging.error("Error: Kv should be an even number.")
        return None

    N = np.zeros([Nv, Nv], dtype=int)
    Kv_half = int(round(Kv/2))
    logging.debug("get_ring_lattice: Kv_half is " + str(Kv_half))

    for i in range(0, Nv):
        for c in range(0, Kv_half):
            j = i + c + 1
            j = j % Nv
            N[i][j] = 1
            N[j][i] = 1

    return N


def get_rewired_lattice(Nv, Kv, pR):
    """
    1. Generate a ring lattice network
    2. Rewire edges

    :param Nv: Number of vertices. int in [0,1000].
    :param pR: Probability of rewiring Kv/2 edges to the right. float in [0,1].
    :param Kv: Average degree. int in [0,Nv-1]; must be an even number.
    :return: N: Matrix as the small-world network
    """
    if type(pR) is not float:
        logging.error("get_rewired_lattice: Error: pR should be a float data type.")
        return None

    if pR < 0 or pR > 1:
        logging.error("get_rewired_lattice: Error: pR should be in [0,1].")
        return None

    N = get_ring_lattice(Nv, Kv)
    logging.debug("get_rewired_lattice: Nv is " + str(Nv))
    logging.debug("get_rewired_lattice: Kv is " + str(Kv))
    logging.debug("get_rewired_lattice: pR is " + str(pR))
    Kv_half = int(round(Kv / 2))
    for i in range(0, Nv):
        for c in range(0, Kv_half):
            j = i + c + 1
            if j >= Nv:
                j = j - Nv
            r1 = rg.random()
            if r1 <= pR:
                N[i][j] = 0
                N[j][i] = 0
                k = None
                while k is None:
                    r2 = rg.random()
                    k = int(round(r2 * (Nv - 1)))
                    if k != i and k != j and N[i][k] == 0:
                        N[i][k] = 1
                        N[k][i] = 1
                    else:
                        k = None
    return N


def find_vertex_at_nearest_distance(DISTANCES, D):
    """
    In case there are multiple vertices at the nearest distance,
    the tie is broken at random, and a single vertex is returned.

    :param DISTANCES: Array containing distances from a source vertex
    :param D: The required distance
    :return:
    """
    v = int(0)  # All vertex IDs are integers
    iv = int(0)  # Index of the vertex v in DISTANCES
    DISTANCES = np.asarray(DISTANCES)
    min_val = (np.abs(DISTANCES - D)).min()
    vertices = np.where(DISTANCES == min_val + D)
    iv = int(np.random.random() * (len(vertices[0]) - 1))
    v = vertices[0][iv]
    return v


def get_fully_connected_network(Nv):
    """
    1. Generate a fully connected network with Nv nodes

    :param Nv: Number of vertices. int in [0,1000].
    :return: N: Matrix as the small-world network
    """
    if type(Nv) is not int:
        logging.error("Error: Nv should be an int data type.")
        return None

    if Nv < 0 or Nv > 1000:
        logging.error("Error: Nv should be in [0,1000].")
        return None

    N = np.ones([Nv, Nv], dtype=int)

    logging.debug("get_fully_connected_network: Nv is " + str(Nv))
    return N


def get_extremist_network(Nx, Dx):
    """
    Consider each edge and with probability Dx wire it.

    :param Nx: Number of extremists
    :param Dx: probability of being a clique
    :return: N: Matrix representing the extremist network
    """
    N = np.zeros([Nx, Nx], dtype=int)

    for i in range(0, Nx):  # Iterate over the upper diagonal only
        for j in range(i+1, Nx):
            r = rg.random()  # r is in [0.0, 1.0)
            if r <= Dx:
                N[i][j] = 1
                N[j][i] = 1
    return N


def combine_M_X1_X2(M, X1, X2):
    """
    Combine the matrices representing the moderates (M),
    positive extremists (X1) and negative extremists (X2),
    into a single matrix, that is the 'world' of all agents.
    :param M: Matrix that is the network of moderates
    :param X1: Matrix that is the network of positive extremists
    :param X2: Matrix that is the network of negative extremists
    :return: W: Matrix that is the world with all agents.
    """
    Nm = M.shape[0]
    Nx1 = X1.shape[0]
    Nx2 = X2.shape[0]
    Nt = Nm + Nx1 + Nx2
    logging.debug("combine_M_X1_X2: The total number of agents in the network are:" + str(Nt))
    W = np.zeros([Nt, Nt], dtype=int)

    for i in range(0, Nm):
        for j in range(0, Nm):
            W[i][j] = M[i][j]

    x_x1 = 0
    y_x1 = 0
    for x in range(Nm, Nm + Nx1):
        for y in range(Nm, Nm + Nx1):
            W[x][y] = X1[x_x1][y_x1]
            y_x1 = y_x1 + 1
        y_x1 = 0
        x_x1 = x_x1 + 1

    x_x2 = 0
    y_x2 = 0
    for x in range(Nm + Nx1, Nt):
        for y in range(Nm + Nx1, Nt):
            W[x][y] = X2[x_x2][y_x2]
            y_x2 = y_x2 + 1
        y_x2 = 0
        x_x2 = x_x2 + 1

    return W


def find_shortest_path_to_other_vertices_using_Matrix_Multiply(N, v0):
    """
    Find the shortest path to vertices in M, from v0.

    :param N: Matrix as a small-world network
    :param v0: Vertex in N from which to find the longest path
    :return: DISTANCES: array containing distances from v0, to other vertices.
    """
    Nv = N.shape[0]
    N__ = np.copy(N)
    DISTANCES = np.full([Nv], -1, dtype=int)
    DISTANCES[v0] = 0
    max_path_len_possible = Nv - 1
    keep_searching = True
    i = 0
    while i <= max_path_len_possible and keep_searching == True:
        N_ = np.copy(N__)
        all_marked = True
        ii = 0
        while ii < Nv:
            if DISTANCES[ii] == -1 and N__[v0][ii] != 0:
                DISTANCES[ii] = i + 1  # shortest path to node ii from v0
            elif DISTANCES[ii] == -1:
                all_marked = False  # not all paths found
            ii = ii + 1
        if all_marked == True:
            keep_searching = False
        N__ = np.dot(N_, N)
        i = i + 1

    return DISTANCES


def connect_XtoM(W, Nm, Nx1, Nx2, Dsep, Sx, pX1, pM1, pX2, pM2):
    """
    Connect the two extremist networks, X1 (X+) and X2 (X-), to the
    Moderate network, M, contained in the World, W. With Nm, total
    moderate agents and, Nx1, total extremist agents in X1, and Nx2,
    total extremist agents in X2.

    Connect X1 and X2 at different starting locations in M, m0x1 and m0x2,
    respectively, that are at a distance controlled by the parameter Desp,
    within a neighbourhood of radius bounded by the parameter radius.

    With probability, pX1 (for X1), and pX2 (for X2), that extremist agents
    attempt to rewire existing connections and the corresponding probabilities,
    pM1 (for X1) and pM2 (for X2), that a moderate agent is available for
    rewiring.

    :param W: The world with all agents.
    :param Nm: Total number of moderate agents in W.
    :param Nx1: Total number of extremist agents in X1.
    :param Nx2: Total number of extremist agents in X2.
    :param Dsep: in [0, 1]. Determines the distance between m+ (m0x1) and m- (m0x2).
    :param Sx: in [0, 1]. Radius of the neighbourhood around that will be explored around m+ (for X1) and m- (for X2).
    :param pX1: in [0, 1]. Probability that an agent in X1 will rewire their connection.
    :param pM1: in [0, 1]. Probability that a moderate agent is available for rewiring with X1.
    :param pX2: in [0, 1]. Probability that an agent in X2 will rewire their connection.
    :param pM2: in [0, 1]. Probability that a moderate agent is available for rewiring with X2.
    :return: W: The world with rewired connections.
    """
    if type(Nm) is not int:
        logging.error("connect_XtoM: Nm should be an int.")
        return None
    if type(Nx1) is not int:
        logging.error("connect_XtoM: Nx1 should be an int.")
        return None
    if type(Nx2) is not int:
        logging.error("connect_XtoM: Nx2 should be an int.")
        return None
    if type(Dsep) is not float:
        logging.error("connect_XtoM: Dsep should be a float.")
        return None
    if type(Sx) is not float:
        logging.error("connect_XtoM: Sx should be a float.")
        return None
    if type(pX1) is not float:
        logging.error("connect_XtoM: pX1 should be a float.")
        return None
    if type(pX2) is not float:
        logging.error("connect_XtoM: pX2 should be a float.")
        return None
    if type(pM1) is not float:
        logging.error("connect_XtoM: pM1 should be a float.")
        return None
    if type(pM2) is not float:
        logging.error("connect_XtoM: pM2 should be a float.")
        return None
    if Nm < 0 or Nm > 1000:
        logging.error("connect_XtoM: Nm should be in [0, 900].")
        return None
    if Nx1 < 0 or Nx1 >= Nm:
        logging.error("connect_XtoM: Nx1 should be in [0, Nm).")
        return None
    if Nx2 < 0 or Nx2 >= Nm:
        logging.error("connect_XtoM: Nx2 should be in [0, Nm).")
        return None
    if W.shape[0] != Nm + Nx1 + Nx2:
        logging.error("connect_XtoM: Total agents in W should be sum of agents in M, X1 and X2.")
        return None
    if Dsep < 0 or Dsep > 1:
        logging.error("connect_XtoM: Dsep should be in [0, 1].")
        return None
    if Sx < 0 or Sx > 1:
        logging.error("connect_XtoM: Sx should be in [0, 1].")
        return None
    if pX1 < 0 or pX1 > 1:
        logging.error("connect_XtoM: pX1 should be in [0, 1].")
        return None
    if pX2 < 0 or pX2 > 1:
        logging.error("connect_XtoM: pX2 should be in [0, 1].")
        return None
    if pM1 < 0 or pM1 > 1:
        logging.error("connect_XtoM: pM1 should be in [0, 1].")
        return None
    if pM2 < 0 or pM2 > 1:
        logging.error("connect_XtoM: pM2 should be in [0, 1].")
        return None

    Nw = W.shape[0]  # Sum count of all agents from the three networks M, X1 and X2
    Km = int(0)  # Calculated below: the average degree of the moderate network
    M = np.zeros([Nm, Nm], dtype=int)  # Agents in M, with all edges, copied from W
    X = np.zeros([Nx1 + Nx2, Nx1 + Nx2], dtype=int)  # combined matrix for agents in X1 and X2, with edges copied only into the upper diagonal, taken from W
    x_type = int(0)  # Value in {1,2} to represent agent in X1 or X2 picked from X
    xy = np.zeros([2], dtype=int)  # edge to be rewired from X1 or X2
    ij = np.zeros([2], dtype=int)  # edge to be rewired from M
    c_xy = int(0)  # random index of an edge in X1 or X2
    rn = int(0)  # ID of a random agent
    i_rn = int(0)  # index into an array of ad-hoc size containing neighbours of a given agent
    m0x1 = int(round(rg.random() * (Nm-1)))  # moderate agent from which agents in X1 start their walk
    m0x2 = int(0)  # moderate agent, at a specified distance from m0x1, from which agents in X2 start their walk
    m = int(0)  # moderate agent from which to take the next step in the walk
    DISTANCES = np.full([Nm], -1, dtype=int)  # Store geodesic from m0x1 to every other node in M.
    Dmax = int(0)  # The largest geodesic found after calculating DISTANCES.
    L = int(0)  # radius of the random walk in M. This value sets the path length of the random walk in M.
    w = int(0)  # Number of steps in the walk, re-evaluated for each xy in X
    Nex = int(0)  # Total undirected edges in X1 and X2
    NEIGHBOURS = None  # used in the random walk to identify a random neighbour

    for i in range(0, Nm):
        for j in range(0, Nm):
            if W[i][j] == 1:
                Km = Km + 1
    Km = int(round(Km / Nm))
    logging.debug("connect_XtoM: Km is " + str(Km))
    L = int(Sx * Nm * Km )  # L is in [0, Nm * (Nm - 1)]
    logging.debug("connect_XtoM: Sx is " + str(Sx))
    logging.debug("connect_XtoM: L is " + str(L))
    # Extract M from W and copy the full matrix
    for i in range(0, Nm):
        for j in range(0, Nm):
            M[i][j] = W[i][j]

    # Extract X1 from W, however, only copy the upper diagonal; and
    # add the total number of edges in X1 to Nex.
    x_W = Nm  # Agents in X1 are placed after agents in M, in W.
    y_W = Nm
    for x in range(0, Nx1):
        for y in range(x, Nx1):
            X[x][y] = W[x_W][y_W]
            y_W = y_W + 1
            if X[x][y] == 1:
                Nex = Nex + 1
        y_W = x_W + 1
        x_W = x_W + 1
    logging.debug("connect_XtoM: completed copying X1 to X and Nex is: " + str(Nex))
    # Now extract X2 from W and copy only the upper diagonal.
    # Add the total number of edges in X2 to Nex.
    x_W = Nm + Nx1  # Agents in X2 are placed after agents in M and X1, in W.
    y_W = Nm + Nx1
    for x in range(Nx1, Nx1 + Nx2):
        for y in range(x, Nx1 + Nx2):
            X[x][y] = W[x_W][y_W]
            y_W = y_W + 1
            if X[x][y] == 1:
                Nex = Nex + 1
        y_W = x_W + 1
        x_W = x_W + 1
    logging.debug("connect_XtoM: completed copying X2 to X and Nex is: " + str(Nex))
    # print(X)
    # Nex is now equal to the the sum of the undirected edges in X1 and X2.
    DISTANCES = find_shortest_path_to_other_vertices_using_Matrix_Multiply(M, m0x1)
    Dmax = DISTANCES.max()
    while Dmax == 0:  # Search to locate a non-isolated node
        m0x1 = int(round(rg.random() * (Nm-1)))
        DISTANCES = find_shortest_path_to_other_vertices_using_Matrix_Multiply(M, m0x1)
        Dmax = DISTANCES.max()
    D = int(round(Dmax * Dsep))
    m0x2 = find_vertex_at_nearest_distance(DISTANCES, D)

    # print(DISTANCES)
    logging.debug("connect_XtoM: Dmax is " + str(Dmax))
    logging.debug("connect_XtoM: m0x1 is " + str(m0x1))
    logging.debug("connect_XtoM: D is " + str(D))
    logging.debug("connect_XtoM: m0x2 is " + str(m0x2))
    logging.debug("connect_XtoM: pX1 is " + str(pX1))
    logging.debug("connect_XtoM: pM1 is " + str(pM1))
    logging.debug("connect_XtoM: px2 is " + str(pX2))
    logging.debug("connect_XtoM: pM2 is " + str(pM2))

    while Nex > 0:  # Consider all edges in X1 and X2, using the combined matrix X
        logging.debug("connect_XtoM: Nex is now: " + str(Nex))
        x_type = 0
        xy = np.zeros([2], dtype=int)
        ij = np.zeros([2], dtype=int)
        c_xy = int(round(rg.random() * (Nex - 1)))  # rg.random() is in [0, 1.0)
        # c_xy is in [0, Nex - 1]
        c = int(0)
        for x in range(0, Nx1 + Nx2):
            for y in range(x + 1, Nx1 + Nx2):
                if X[x][y] == 1 and c < c_xy:
                    c = c + 1
                elif X[x][y] == 1 and c == c_xy:
                    c = c + 1  # since the for loop keeps running
                    if x < Nx1:
                        x_type = 1
                    else:
                        x_type = 2
                    xy[0] = x
                    xy[1] = y
                    logging.debug("connect_XtoM: x,y is: " + str(xy[0]) + ", " + str(xy[1]))
                    X[x][y] = 0
                    Nex = Nex - 1
        r2 = rg.random()
        if (x_type == 1 and r2 <= pX1) or (x_type == 2 and r2 <= pX2):  # consider rewiring xy
            if (x_type == 1 and rg.random() <= pM1) or (x_type == 2 and rg.random() <= pM2):  # Q: Evaluate pM before starting the walk. Is this valid?
                c = int(0)
                r4 = rg.random()
                w = int(round(r4 * (L - 1)))  # w: the number of steps in the walk will be in [0, L - 1]
                if x_type == 1:
                    m = m0x1
                elif x_type == 2:
                    m = m0x2

                while w > 0:  # continue walk for w steps
                    NEIGHBOURS = None
                    NEIGHBOURS = np.where(W[m] == 1)
                    i_rn = int(round(rg.random() * (len(NEIGHBOURS[0]) - 1)))
                    rn = NEIGHBOURS[0][i_rn]
                    w = w - 1
                    m = rn

                # next: find a random neighbour of m and store edge to that agent as ij
                ij[0] = m
                NEIGHBOURS = None
                NEIGHBOURS = np.where(W[m] == 1)
                i_rn = int(round(rg.random() * (len(NEIGHBOURS[0]) - 1)))
                ij[1] = NEIGHBOURS[0][i_rn]
                # print(NEIGHBOURS[0])
                logging.debug("connect_XtoM: i,j is: " + str(ij[0]) + ", " + str(ij[1]))
                # print("before rewiring, W, is: ")
                # print(W)

                if rg.random() <= 0.5:  # remove edge xy
                    logging.debug("removing x,y")
                    W[xy[0] + Nm][xy[1] + Nm] = 0
                    W[xy[1] + Nm][xy[0] + Nm] = 0
                else:  # remove edge ij
                    logging.debug("removing i,j")
                    W[ij[0]][ij[1]] = 0
                    W[ij[1]][ij[0]] = 0

                if rg.random() <= 0.5:  # add edge xi
                    logging.debug("add x,i")
                    W[xy[0] + Nm][ij[0]] = 1
                    W[ij[0]][xy[0] + Nm] = 1
                else:  # add edge yi
                    logging.debug("add y,i")
                    W[xy[1] + Nm][ij[0]] = 1
                    W[ij[0]][xy[1] + Nm] = 1

                # print("after rewiring, W, is: ")
                # print(W)
                # print("\n\n")

    return W


def generate_network(Nm, Km, pR, Nx1, Nx2):
    """

    :param Nm: Total number of moderate agents
    :param Km: Average degree of the moderate network
    :param pR: probability with which to re-wire the initial ring-lattice moderate network
    :param pe: the proportion of extremists in the population
    :param Nx1: Total number of positive extremists
    :param Nx2: Total number of negative extremists
    :return: The World, W, as an adjacency matrix with all agents.
    """
    # Social Network Parameters:
    #  Nw: Total agents in the network including extremists
    # Parameters for moderates
    #  Nm = int(0) # Total number of moderates agents
    #  Km: Average degree of small-world network of moderates
    #  pR: rewiring probability of initial ring-lattice network of moderates
    # Parameters for positive extremists
    tot_extremists = int(0)
    Dx1 = float(1) # density of edges between positive extremists normalized in the range [0, 1]
    pX1 = float(1) # proba
    # bility that a positive extremist will rewire with a moderate
    pM1 = float(1) # probability that a moderate will rewire with a positive extremist
    # Parameters for negative extremists
    Dx2 = float(1) # density of edges between negative extremists normalized in the range [0, 1]
    pX2 = float(1) # probability that a negative extremist will rewire with a moderate
    pM2 = float(1) # probability that a moderate will rewire with a negative extremist
    # Parameters for connecting extremists with moderates
    Dsep = float(1)  # Dsep in [0, 1]. Determines the distance between first point of contacts in the moderate network of each of the extremist community
    Sx = float(0.5)  # Sx in [0, 1]  Diffusion of each extremist community in the moderate network

    # Create the social network based on the value of pe:
    #  Create a network of moderates:
    logging.debug(("generate_network: Nx1 is " + str(Nx1)))
    logging.debug(("generate_network: Nx2 is " + str(Nx2)))

    M = get_rewired_lattice(Nm, Km, pR)
    logging.debug("generate_network: Created moderate network M.")
    #  Create extremist networks
    X1 = get_extremist_network(Nx1, Dx1)
    logging.debug("generate_network: X1 created with Dx1="+str(Dx1))
    # print(X1)
    X2 = get_extremist_network(Nx2, Dx2)
    logging.debug("generate_network: Created extremist network X2 created with Dx2="+str(Dx2))
    # print(X2)
    # Connect the networks:
    # print(W)
    W = combine_M_X1_X2(M, X1, X2) # create a composite matrix
    logging.debug("generate_network: Combined X1, X2 and M")
    W = connect_XtoM(W, Nm, Nx1, Nx2, Dsep, Sx, pX1, pM1, pX2, pM2)  # Add edges between extremist and moderate agents
    logging.debug("generate_network: Connected X1, X2 and M")
    # print("After rewiring W is:")
    # print(W)
    # We now have a world, W, with agents that are all connected based on the input parameters, and now,
    return W


def main():
    logging.basicConfig(level=logging.ERROR)
    topology = str(sys.argv[1])
    #topology = "ring-lattice"
    pe = float(sys.argv[2]) # proportion of extremists in the population
    #pe = float(0.3)
    W = None  # The adjacency matrix that represents the world with all agents
    pR = float(0) # Re-wiring probability of the inital ring-lattice moderate network
    Nm = int(1000) # Total number of moderate agents
    Nx1 = int(0)  # Total number of positive extremists
    Nx2 = int(0)  # Total number of negative extremists
    Km = int(800) # Average degree of the moderate network
    tot_extremists = int(Nm*pe) # Total number of extremists in the population
    if is_even(tot_extremists) == False:
        tot_extremists += 1
    logging.debug("main: tot_extremists is: " + str(tot_extremists))
    Nx1 = int(tot_extremists/2)
    Nx2 = Nx1
    if topology == "ring-lattice":
        pR = float(0)
    elif topology == "small-world":
        pR = float(0.03)
    elif topology == "random":
        pR = float(1)
    export_filename = topology+"_network_"+str(pe)+"_K400.csv"

    # Create target Directory if don't exist
    network_directory = "./networks/"
    if not os.path.exists(network_directory):
        os.mkdir(network_directory)

    start = timer()
    W = generate_network(Nm, Km, pR, Nx1, Nx2)
    end = timer()
    total_time = end - start
    print("Total time to generate the network was %.2f" % total_time+" seconds.\n\n")
    export_world_as_CSV(W, Nm, Nx1, Nx2, export_filename)
    return


main()
sys.exit(0)
