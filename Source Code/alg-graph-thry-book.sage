"""
Functions to deal with oriented graphs as in Biggs's and Baker's works, for
use in [2]. 
Added in 2016: functions to construct expander graphs from linear error-correcting block codes. Also, functions to compute combinatorial calculus functions.
Added in 2017: combinatorial calculus functions for simplicies.

REFERENCES:
    [0] M. Baker and S. Norine, "Riemann-Roch and Abel-Jacobi
    theory on a finite graph", available:
    https://arxiv.org/abs/math/0608360
    [1] N. Biggs, "Chip firing and the critical group of a 
    graph," J. Algeb. Comb. 9 (1999)25-45.
    [2] D. Joyner and C. Melles, "Adventures in graph theory,"
    Birkhauser, to appear.

last modified 2017-07-17

Copyright 2014-2017 - David Joyner wdjoyner@gmail.com

Licensed modified BSD

"""

def tree_component_plus(Gamma, e, T, eo):
    """
    This computes the positive component T_e^+ of a 
    spanning tree T in a graph Gamma, with e in T.

    INPUT:
        Gamma - graph
        e  - edge of Gamma
        T  - spanning tree of Gamma  
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: T = Gamma.min_spanning_tree()
        sage: E = Gamma.edges()
        sage: e = E[5]
        sage: e in T
        True
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: tree_component_plus(Gamma, e, T, eo)
        Subgraph of (): Graph on 1 vertex
        sage: tree_component_plus(Gamma, e, T, eo).vertices()
        [a + 2]
        sage: tree_component_minus(Gamma, e, T, eo)
        Subgraph of (): Graph on 8 vertices

    """
    E = Gamma.edges()
    GammaT = Graph(T)
    if not(e in T):
        return T
    k = E.index(e)
    GammaT.delete_edge(E[k])
    GammaTcc = GammaT.connected_components_subgraphs()
    #print GammaTcc
    GammaT1 =  GammaTcc[0]
    GammaT2 =  GammaTcc[1]
    V1 = GammaT1.vertices()
    V2 = GammaT2.vertices()
    if eo[k] == 1:
        v = e[0]
    else:
        v = e[1]
    if v in V1:
        return GammaT1
    else:
        return GammaT2

def tree_component_minus(Gamma, e, T, eo):
    """
    This computes the positive component T_e^+ of a 
    spanning tree T in a graph Gamma, with e in T.

    INPUT:
        Gamma - graph
        e  - edge of Gamma
        T  - spanning tree of Gamma  
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: T = Gamma.min_spanning_tree()
        sage: E = Gamma.edges()
        sage: e = E[5]
        sage: e in T
        True
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: tree_component_plus(Gamma, e, T, eo)
        Subgraph of (): Graph on 1 vertex
        sage: tree_component_plus(Gamma, e, T, eo).vertices()
        [a + 2]
        sage: tree_component_minus(Gamma, e, T, eo).vertices()
        [0, 1, 2, a, a + 1, 2*a, 2*a + 1, 2*a + 2]

    """
    E = Gamma.edges()
    GammaT = Graph(T)
    if not(e in T):
        return T
    k = E.index(e)
    GammaT.delete_edge(E[k])
    GammaTcc = GammaT.connected_components_subgraphs()
    #print GammaTcc
    GammaT1 =  GammaTcc[0]
    GammaT2 =  GammaTcc[1]
    V1 = GammaT1.vertices()
    V2 = GammaT2.vertices()
    if eo[k] == -1:
        v = e[0]
    else:
        v = e[1]
    if v in V1:
        return GammaT1
    else:
        return GammaT2

def incidence_value(Gamma, v, e, eo):
    """
    This computes the incidence value of a vertex and edge of 
    a graph Gamma with edge-orientation vector eo. 

    INPUT:
        Gamma - graph
        v  - vertex of Gamma
        e  - edge of Gamma  
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
 
    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: E = Gamma.edges()
        sage: V = Gamma.vertices()
        sage: eo = [1]*len(E)
        sage: incidence_value(Gamma, V[2], E[3], eo)
        0
        sage: incidence_value(Gamma, V[8], E[3], eo)
        -1

    """
    E = Gamma.edges()
    if v in e:
        if v == e[0]:
            k = E.index(e)
            return eo[k]
        elif v == e[1]:
            k = E.index(e)
            return -eo[k]
        else:
            return 0
    return 0

def incidence_matrix(Gamma, eo):
    """
    This computes the incidence matrix (whose rows are indexed by edges
    and whose columns are indexed by vertices) of a graph Gamma with 
    edge-orientation vector eo. The ordering of the edges and of the vertices is the same 
    as Sage's vertices and edges methods.

    INPUT:
        Gamma - graph
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: E = Gamma.edges()
        sage: V = Gamma.vertices()
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: incidence_matrix(Gamma, eo)
        [ 1 -1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
        [-1  0  0  0 -1 -1  1  0  0  0  0  0  0  0  0  0  0  0]
        [ 0  1  0  0  1  0  0 -1  1  0  0  0  0  0  0  0  0  0]
        [ 0  0  0  0  0  0  0  1  0  1 -1 -1  0  0  0  0  0  0]
        [ 0  0 -1  0  0  0  0  0  0 -1  0  0  1 -1  0  0  0  0]
        [ 0  0  0  0  0  1  0  0  0  0  1  0 -1  0  1  0  0  0]
        [ 0  0  0  0  0  0 -1  0  0  0  0  0  0  0 -1  1 -1  0]
        [ 0  0  0  0  0  0  0  0 -1  0  0  1  0  0  0 -1  0 -1]
        [ 0  0  0 -1  0  0  0  0  0  0  0  0  0  1  0  0  1  1]

        sage: IM = incidence_matrix(Gamma, eo)
        sage: IM1 = IM.delete_rows([0])
        sage: Proj1 = IM1.transpose()*(IM1*IM1.transpose())^(-1)*IM1
        sage: Proj1.image() # row space
        Vector space of degree 18 and dimension 8 over Rational Field
        Basis matrix:
        [ 1  0  0  0  1  0  0  0  0  0 -1  0  1  0  0 -1  1  0]
        [ 0  1  0  0  1  0  0  0  0  1 -1  0  0  0  0 -1  0 -1]
        [ 0  0  1  0  0  0  0  0  0  1  0  0 -1  1  0  0  0  0]
        [ 0  0  0  1  0  0  0  0  0  0  0  0  0 -1  0  0 -1 -1]
        [ 0  0  0  0  0  1  0  0  0  0  1  0 -1  0  1  0  0  0]
        [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  1 -1  1  0]
        [ 0  0  0  0  0  0  0  1  0  1 -1 -1  0  0  0  0  0  0]
        [ 0  0  0  0  0  0  0  0  1  0  0 -1  0  0  0  1  0  1]
        sage: cut_space(Gamma, eo)
        Vector space of degree 18 and dimension 8 over Rational Field
        Basis matrix:
        [ 1  0  0  0  1  0  0  0  0  0 -1  0  1  0  0 -1  1  0]
        [ 0  1  0  0  1  0  0  0  0  1 -1  0  0  0  0 -1  0 -1]
        [ 0  0  1  0  0  0  0  0  0  1  0  0 -1  1  0  0  0  0]
        [ 0  0  0  1  0  0  0  0  0  0  0  0  0 -1  0  0 -1 -1]
        [ 0  0  0  0  0  1  0  0  0  0  1  0 -1  0  1  0  0  0]
        [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  1 -1  1  0]
        [ 0  0  0  0  0  0  0  1  0  1 -1 -1  0  0  0  0  0  0]
        [ 0  0  0  0  0  0  0  0  1  0  0 -1  0  0  0  1  0  1]
        sage: B = incidence_matrix(Gamma, eo) 
        sage: B*B.transpose() == Gamma.laplacian_matrix()
        True

        sage: V = range(10)
        sage: E = [(0,1),(0,2),(0,5),(1,2),(1,3),(1,5),(2,4),(2,5),(3,4),(3,5),(3,6),\
        ....: (4,5),(4,7),(5,6),(5,7),(5,8),(6,7),(6,8),(7,8),(8,9)]
        sage: Gamma = Graph([V,E]) # Seifrot graph of Kaballah mysticism
        sage: eo = 20*[1]
        sage: incidence_matrix(Gamma, eo)
        [ 1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
        [-1  0  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
        [ 0 -1  0 -1  0  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0]
        [ 0  0  0  0 -1  0  0  0  1  1  1  0  0  0  0  0  0  0  0  0]
        [ 0  0  0  0  0  0 -1  0 -1  0  0  1  1  0  0  0  0  0  0  0]
        [ 0  0 -1  0  0 -1  0 -1  0 -1  0 -1  0  1  1  1  0  0  0  0]
        [ 0  0  0  0  0  0  0  0  0  0 -1  0  0 -1  0  0  1  1  0  0]
        [ 0  0  0  0  0  0  0  0  0  0  0  0 -1  0 -1  0 -1  0  1  0]
        [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0 -1 -1  1]
        [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1]

    """
    E = Gamma.edges()
    V = Gamma.vertices()
    IG = [[incidence_value(Gamma, v, e, eo) for v in V] for e in E]
    #print IG
    return matrix(QQ, IG).transpose()

def cycle_space(Gamma, eo, F = QQ):
    """
    This computes the cycle space Z of Gamma.

    INPUT:
        Gamma - graph
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: E = Gamma.edges()
        sage: V = Gamma.vertices()
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: cycle_space(Gamma, eo)
        Vector space of degree 18 and dimension 10 over Rational Field
        Basis matrix:
        [ 1  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0 -1  0]
        [ 0  1  0  1  0  0  0  0 -1  0  0  0  0  0  0  0  0  1]
        [ 0  0  1 -1  0  0  0  0  0  0  0  0  0 -1  0  0  0  0]
        [ 0  0  0  0  1  0  1  0 -1  0  0  0  0  0  0  0 -1  1]
        [ 0  0  0  0  0  1  1  0  0  0  0  0  0  0 -1  0  0  0]
        [ 0  0  0  0  0  0  0  1  1  0  0  1  0  0  0  0  0  0]
        [ 0  0  0  0  0  0  0  0  0  1  0  1  0 -1  0  0  0  1]
        [ 0  0  0  0  0  0  0  0  0  0  1 -1  0  0 -1  0  1 -1]
        [ 0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  0 -1  0]
        [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1 -1]
    """
    m = len(Gamma.edges())
    n = len(Gamma.vertices())
    MS = MatrixSpace(F, m, n)
    IG = MS(incidence_matrix(Gamma, eo).transpose())
    return IG.kernel()

def cocycle_space(Gamma, eo):
    return cut_space(Gamma, eo)

def cut_space(Gamma, eo):
    """
    This computes the cut space B of Gamma.

    INPUT:
        Gamma - graph
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: E = Gamma.edges()
        sage: V = Gamma.vertices()
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: cut_space(Gamma, eo)
        Vector space of degree 18 and dimension 8 over Rational Field
        Basis matrix:
        [ 1  0  0  0  1  0  0  0  0  0 -1  0  1  0  0 -1  1  0]
        [ 0  1  0  0  1  0  0  0  0  1 -1  0  0  0  0 -1  0 -1]
        [ 0  0  1  0  0  0  0  0  0  1  0  0 -1  1  0  0  0  0]
        [ 0  0  0  1  0  0  0  0  0  0  0  0  0 -1  0  0 -1 -1]
        [ 0  0  0  0  0  1  0  0  0  0  1  0 -1  0  1  0  0  0]
        [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  1 -1  1  0]
        [ 0  0  0  0  0  0  0  1  0  1 -1 -1  0  0  0  0  0  0]
        [ 0  0  0  0  0  0  0  0  1  0  0 -1  0  0  0  1  0  1]
        sage: v = cut_space(Gamma, eo).random_element()
        sage: w = cycle_space(Gamma, eo).random_element()
        sage: v.inner_product(w)
        0
    """
    Z = cycle_space(Gamma, eo)
    return Z.complement()

def incidence_matrix_tree_inverse_value(Gamma, eo, x, T, e, v):
    """
    This computes the matrix entries of the inverse of the truncated 
    incidence matrix, D(x,T)^(-1), with edge-orientation vector eo, 
    associated to a vertex x in a spanning tree T.

    INPUT:
        Gamma - graph
        x,v  - distinct vertices of Gamma
        T  - spanning tree of Gamma  
        e  - an edge of Gamma in T
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: E = Gamma.edges()
        sage: V = Gamma.vertices()
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: T = Gamma.min_spanning_tree()
        sage: x = V[2]
        sage: [incidence_matrix_tree_inverse_value(Gamma, eo, x, T, e, v) for v in V if v<>x]
        [0, 0, 0, 0, 1, 0, 0, 0]
        sage: e = E[0]
        sage: e in T
        True
        sage: [incidence_matrix_tree_inverse_value(Gamma, eo, x, T, e, v) for v in V if v<>x]
        [0, -1, 0, 0, -1, -1, 0, 0]

    This is consistent with the inverse matrix example for incidence_matrix_tree(Gamma, eo, x, T).

    """
    VT_e_plus = tree_component_plus(Gamma, e, T, eo).vertices()
    VT_e_minus = tree_component_minus(Gamma, e, T, eo).vertices()
    if v in VT_e_plus and x in VT_e_minus:
        return 1
    if x in VT_e_plus and v in VT_e_minus:
        return -1
    return 0

def incidence_matrix_tree(Gamma, eo, x, T):
    """
    This computes the truncated incidence matrix D(x,T) with 
    edge-orientation vector eo, associated to a vertex x
    in a spanning tree T.

    INPUT:
        Gamma - graph
        x  - vertex of Gamma
        T  - spanning tree of Gamma  
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: E = Gamma.edges()
        sage: V = Gamma.vertices()
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: T = Gamma.min_spanning_tree()
        sage: v = V[2]
        sage: incidence_matrix_tree(Gamma, eo, v, T)
        [ 1 -1  0  0  0  0  0  0]
        [-1  0  0  0  0  0  0  0]
        [ 1  0  0 -1  0  0  0  0]
        [ 1  0  0  0  0  0  0 -1]
        [ 0 -1  0  0  1  0  0  0]
        [ 0  1  0  0  0 -1  0  0]
        [ 0  0  1  0  0  0  0  0]
        [ 0  0  0  0  0  0 -1  0]
        sage: DvT^(-1)
        [ 0 -1  0  0  0  0  0  0]
        [-1 -1  0  0  0  0  0  0]
        [ 0  0  0  0  0  0  1  0]
        [ 0 -1 -1  0  0  0  0  0]
        [-1 -1  0  0  1  0  0  0]
        [-1 -1  0  0  0 -1  0  0]
        [ 0  0  0  0  0  0  0 -1]
        [ 0 -1  0 -1  0  0  0  0]

    """
    L = []
    E = Gamma.edges()
    V = Gamma.vertices()
    IM = incidence_matrix(Gamma, eo).transpose()
    i = V.index(x)
    DvT = IM.delete_columns([i])
    EnotT = [e for e in E if not(e in T)]
    for e in EnotT:
        L.append(E.index(e))
    DvT = DvT.delete_rows(L)
    return DvT

def fund_vertex_cut_char_fcn(Gamma, v, e, eo):
    """
    This computes the characteristic function of the fundamental 
    cut defined by a vertex v, evaluated at an edge e of a graph Gamma with 
    edge-orientation vector eo. 

    INPUT:
        Gamma - graph
        v  - vertex of Gamma
        e  - edge of Gamma  
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    OUTPUT:
        1, if v is the head of e (defined by eo)
        -1, if v is the tail of e (defined by eo)
        0, otherwise.

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: E = Gamma.edges()
        sage: V = Gamma.vertices()
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: v = V[2]
        sage: e = E[1]
        sage: fund_vertex_cut_char_fcn(Gamma, v, e, eo)
        1
        sage: [fund_vertex_cut_char_fcn(Gamma, V[2], e, eo) for e in E]
        [0, 1, 0, 0, 1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        # 3rd column of incidence_matrix(Gamma, eo)
        sage: [fund_vertex_cut_char_fcn(Gamma, v, E[1], eo) for v in V]
        [-1, 0, 1, 0, 0, 0, 0, 0, 0]
        # 2nd row of incidence_matrix(Gamma, eo)

    """
    E = Gamma.edges()
    V = Gamma.vertices()
    k = E.index(e)
    IG = incidence_matrix(Gamma, eo).transpose()
    char_fcn = IG.rows()[k]
    j = V.index(v)
    return char_fcn[j]

def fund_cut_char_fcn(Gamma, v, e, eo):
    """
    This computes the characteristic function of the fundamental 
    cut defined by a subset U of vertices, evaluated at an edge e of a graph Gamma with 
    edge-orientation vector eo. The cut defined by U is the set of edges having exactly 
    one vertex in U.

    INPUT:
        Gamma - graph
        U  - a list of vertices of Gamma
        e  - edge of Gamma  
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    OUTPUT:
        1, if v is the head of e (defined by eo)
        -1, if v is the tail of e (defined by eo)
        0, otherwise.

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: E = Gamma.edges()
        sage: V = Gamma.vertices()
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: U = [V[2],V[4]]
        sage: e = E[1]
        sage: fund_cut_char_fcn(Gamma, U, e, eo)
        1
        sage: [fund_cut_char_fcn(Gamma, U, e, eo) for e in E]
        [0, 1, -1, 0, 1, 0, 0, -1, 1, -1, 0, 0, 1, -1, 0, 0, 0, 0]

    """
    L = [fund_vertex_cut_char_fcn(Gamma, v, e, eo) for v in U]
    return sum(L)

def fund_cycle_char_fcn(Gamma, q, e, eo):
    """
    This computes the characteristic function of the fundamental 
    cycle q, evaluated at an edge e of a graph Gamma with 
    edge-orientation vector eo. 

    INPUT:
        Gamma - graph
        q  - fundamental cycle of Gamma
        e  - edge of Gamma  
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    OUTPUT:
        1, if v is the head of e (defined by eo)
        -1, if v is the tail of e (defined by eo)
        0, otherwise.

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: Q = Gamma.cycle_basis()
        sage: V = Gamma.vertices()
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: e = E[1]
        sage: [fund_cycle_char_fcn(Gamma, q, e, eo) for q in Q]
        [0, 0, 0, 0, 0, 0, -1, 0, -1, -1]

    """
    Q = Gamma.cycle_basis()
    Qset = [Set(qq) for qq in Q]
    if Set(q) in Qset:
        E = Gamma.edges()
        V = Gamma.vertices()
        sgn = eo[E.index(e)]
        if e[0] in q and e[1] in q:
            i = q.index(e[0])
            j = q.index(e[1])
            lq = len(q)
            #print e, i, j, lq
            if i-j==1*sgn or i-j==(lq-1)*sgn:
                return 1
            if i-j==(-1)*sgn or i-j==-(lq-1)*sgn:
                return -1 
    return 0

def fund_cycle(Gamma, T, e, eo):
    """
    This computes the fundamental cycle C(T,F) associated to a tree T
    and an edge e of a graph Gamma, not in T, with edge-orientation  
    vector eo. 

    INPUT:
        Gamma - graph
        T  - a tree of Gamma
        e  - edge of Gamma not in T
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    OUTPUT:
        C(T,f)

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: Q = Gamma.cycle_basis()
        sage: V = Gamma.vertices()
        sage: E = Gamma.edges()
        sage: eo = [1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1]
        sage: e = E[9]
        sage: T = Gamma.min_spanning_tree()
        sage: e in T
        False
        sage: fund_cycle(Gamma, T, e, eo)
        [(a, 2, None), (2, 0, None), (0, a + 1, None), (a + 1, a, None)]

    """
    V = Gamma.vertices()
    i = V.index(e[0])
    j = V.index(e[1])
    P1 = Graph(T).shortest_path(V[i],V[j])
    N = len(L)
    P2 = [(L[i],L[i+1],None) for i in range(N-1)]+[(L[N-1],L[0],None)]
    return P2

def laplacian_reduced(Gamma, eo, s):
    """
    This computes the reduced (vertex) Laplacian of a graph Gamma with 
    edge-orientation vector eo, where the row column associated to s has been removed.

    INPUT:
        Gamma - graph having N vertices
        s - a (source) vertex of Gamma
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    OUTPUT:
        (N-1)x(N-1) matrix of integers which agrees with the usual adjacency 
        matrix Laplacian except for the entries associated to the vertex s.

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: eo = [1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1]
        sage: s = V[2]
        sage: LR = laplacian_reduced(Gamma, eo, s)
        sage: LR
        [ 4 -1  0 -1  0  0  0 -1]
        [-1  4  0  0 -1 -1  0  0]
        [ 0  0  4 -1 -1  0 -1  0]
        [-1  0 -1  4 -1  0  0 -1]
        [ 0 -1 -1 -1  4 -1  0  0]
        [ 0 -1  0  0 -1  4 -1 -1]
        [ 0  0 -1  0  0 -1  4 -1]
        [-1  0  0 -1  0 -1 -1  4]
        sage: LR.dimensions()
        (8, 8)
        sage: A = ZZ^8/LR.image()
        sage: A.cardinality()
        11664
        sage: Gamma.spanning_trees_count()
        11664
        sage: det(LR)
        11664
        sage: factor(11664)
        2^4 * 3^6

    """
    L = Gamma.laplacian_matrix()
    V = Gamma.vertices()
    i = V.index(s)
    L1 = L.delete_rows([i])
    L2 = L1.delete_columns([i])
    return matrix(ZZ, L2)

def laplacian_reduced_punctured(Gamma, eo, s, p, q):
    """
    This computes the reduced, punctured (vertex) Laplacian of a graph Gamma with 
    edge-orientation vector eo, where row for p, s are removed and columns for q,s
    are removed. It's determinant is equal to (pq||s), in Biggs notation [B], section 13-14.


    INPUT:
        Gamma - graph having N vertices
        s - a (source) vertex of Gamma
        p.q - two vertices, each distinct from s
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    OUTPUT:
        (N-2)x(N-2) matrix of integers which agrees with the usual
        adjacency matrix Laplacian except for the entries associated to the
        vertex s.

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: eo = [1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1]
        sage: V = Gamma.vertices()
        sage: s = V[2]; p = V[4]; q = V[6]
        sage: LRP = laplacian_reduced_punctured(Gamma, eo, s, p, q)
        sage: LRP
        [ 4 -1  0 -1  0  0 -1]
        [-1  4  0  0 -1 -1  0]
        [ 0  0  4 -1 -1  0  0]
        [-1  0 -1  4 -1  0 -1]
        [ 0 -1  0  0 -1  4 -1]
        [ 0  0 -1  0  0 -1 -1]
        [-1  0  0 -1  0 -1  4]
        sage: LRP.dimensions()
        (7, 7)
        sage: det(LRP)
        2592

    REFERENCES:
        Biggs, Norman, "Algebraic potential theory on graphs", BLMS, 1997.

    """
    LR = laplacian_reduced(Gamma, eo, s)
    V = Gamma.vertices()
    i = V.index(p)
    j = V.index(q)
    L1 = LR.delete_rows([i])
    L2 = L1.delete_columns([j])
    return matrix(ZZ, L2)
    
def potential_induced(Gamma, e0, j, s, p, q):
    """
    Returns the potential function phi:V -> QQ which is induced by a source
    c having magnitude j, with input p and output q and satisfies phi(q)=0.
    Such a phi satisfies the "differential equation" 
                            D^t*phi = Pc,
    where 
      * D is the incidence matrix, and
      * P is the orthogonal projection from the space, C^1(Gamma, RR), 
        of functions on E (the edges of Gamma) to the cut space B.
    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: e0 = [1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1]
        sage: V = Gamma.vertices()
        sage: s = V[2]; p = V[4]; q = V[6]
        sage: j = 5
        sage: potential_induced(Gamma, e0, j, s, p, q)
        2/45

    """
    LR = laplacian_reduced(Gamma, eo, s)
    LRP = laplacian_reduced_punctured(Gamma, eo, s, p, q)
    return det(LRP)/(det(LR)*j)

def number_of_2_trees(Gamma, e0, v, w):
    """
    Returns the number of 2-trees for which one component contains
    p and the other contains q, (p||q).

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: eo = [1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1]
        sage: V = Gamma.vertices()
        sage: v = V[4]; w = V[6]
        sage: number_of_2_trees(Gamma, e0, v, w)
        6480

    """
    LRP = laplacian_reduced_punctured(Gamma, eo, w, v, v)
    return det(LRP)

def critical_group(Gamma):
    """
    Returns the co-kernel of a reduced Laplacian
    of the graph Gamma.

    EXAMPLES:
        sage: WG6 = graphs.WheelGraph(6)
        sage: critical_group(WG6)
        Finitely generated module V/W over Integer Ring with invariants (11, 11)
        sage: G = graphs.PaleyGraph(9); G
        Paley graph with parameter 9: Graph on 9 vertices
        sage: critical_group(G)
        Finitely generated module V/W over Integer Ring with invariants (6, 6, 18, 18)
        sage: 18*18*6*6
        11664
        sage: G = graphs.CompleteGraph(4); G
        Complete graph: Graph on 4 vertices
        sage: critical_group(G)
        Finitely generated module V/W over Integer Ring with invariants (4, 4)
        sage: WG4 = graphs.WheelGraph(4)
        sage: WG4.delete_edge((2,3))
        sage: WG4.add_vertex(5)
        sage: WG4.add_edges([(2,5),(3,5)])
        sage: critical_group(WG4)
        Finitely generated module V/W over Integer Ring with invariants (24)
        sage: WG4 = graphs.WheelGraph(4)
        sage: critical_group(WG4)
        Finitely generated module V/W over Integer Ring with invariants (4, 4)

        sage: G = DihedralGroup(4)
        sage: G.gens()
        [(1,2,3,4), (1,4)(2,3)]
        sage: S = G.gens()+[G.gens()[0]^(-1)]; S
        [(1,2,3,4), (1,4)(2,3), (1,4,3,2)]
        sage: A = G.cayley_graph(generators=S,simple=True).adjacency_matrix()
        sage: Gamma = Graph(A)
        sage: critical_group(Gamma)
        Finitely generated module V/W over Integer Ring with invariants (2, 8, 24)

    """
    L = Gamma.laplacian_matrix()
    L1 = L.delete_rows([0])
    L2 = L1.delete_columns([0])
    n = L2.dimensions()[0]
    return ZZ^n/L2.image()

def laplacian_vertex(Gamma, eo):
    """
    This computes the vertex Laplacian of a graph Gamma with 
    edge-orientation vector eo. 

    INPUT:
        Gamma - graph having N vertices
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    OUTPUT:
        NxN matrix of 1, -1, 0's (except on the diagonal) which agrees 
        with the usual adjacency matrix Laplacian 

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: eo = [1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1]
        sage: laplacian_vertex(Gamma, eo)
        [ 4 -1 -1  0 -1  0  0  0 -1]
        [-1  4 -1  0  0 -1 -1  0  0]
        [-1 -1  4 -1  0  0  0 -1  0]
        [ 0  0 -1  4 -1 -1  0 -1  0]
        [-1  0  0 -1  4 -1  0  0 -1]
        [ 0 -1  0 -1 -1  4 -1  0  0]
        [ 0 -1  0  0  0 -1  4 -1 -1]
        [ 0  0 -1 -1  0  0 -1  4 -1]
        [-1  0  0  0 -1  0 -1 -1  4]
        sage: Gamma.laplacian_matrix()
        [ 4 -1 -1  0 -1  0  0  0 -1]
        [-1  4 -1  0  0 -1 -1  0  0]
        [-1 -1  4 -1  0  0  0 -1  0]
        [ 0  0 -1  4 -1 -1  0 -1  0]
        [-1  0  0 -1  4 -1  0  0 -1]
        [ 0 -1  0 -1 -1  4 -1  0  0]
        [ 0 -1  0  0  0 -1  4 -1 -1]
        [ 0  0 -1 -1  0  0 -1  4 -1]
        [-1  0  0  0 -1  0 -1 -1  4]
        sage: laplacian_vertex(Gamma, eo).kernel()
        Vector space of degree 9 and dimension 1 over Rational Field
        Basis matrix:
        [1 1 1 1 1 1 1 1 1]


    """
    D = incidence_matrix(Gamma, eo)
    return D.transpose()*D
    
def laplacian_edge(Gamma, eo):
    """
    This computes the edge Laplacian of a graph Gamma with 
    edge-orientation vector eo. 

    INPUT:
        Gamma - graph having N vertices
        eo - a vector of 1's and -1's whose length is the number of edges in Gamma
             (ie, the size of Gamma, M)

    OUTPUT:
        MxM matrix of 1, -1, 0's which agrees with the edge
        adjacency matrix Laplacian 

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(9)
        sage: eo = [1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1]
        sage: laplacian_edge(Gamma, eo)
        [ 2  1 -1  1 -1  1  1  0  0  0  0  0  0  0  0  0  0  0]
        [ 1  2 -1  1  1  0  0  1 -1  0  0  0  0  0  0  0  0  0]
        [-1 -1  2 -1  0  0  0  0  0 -1  0  0  1 -1  0  0  0  0]
        [ 1  1 -1  2  0  0  0  0  0  0  0  0  0 -1  0  0  1 -1]
        [-1  1  0  0  2 -1 -1  1 -1  0  0  0  0  0  0  0  0  0]
        [ 1  0  0  0 -1  2  1  0  0  0  1  0 -1  0 -1  0  0  0]
        [ 1  0  0  0 -1  1  2  0  0  0  0  0  0  0  1 -1  1  0]
        [ 0  1  0  0  1  0  0  2 -1  1 -1 -1  0  0  0  0  0  0]
        [ 0 -1  0  0 -1  0  0 -1  2  0  0 -1  0  0  0 -1  0  1]
        [ 0  0 -1  0  0  0  0  1  0  2 -1 -1 -1  1  0  0  0  0]
        [ 0  0  0  0  0  1  0 -1  0 -1  2  1 -1  0 -1  0  0  0]
        [ 0  0  0  0  0  0  0 -1 -1 -1  1  2  0  0  0  1  0 -1]
        [ 0  0  1  0  0 -1  0  0  0 -1 -1  0  2 -1  1  0  0  0]
        [ 0  0 -1 -1  0  0  0  0  0  1  0  0 -1  2  0  0 -1  1]
        [ 0  0  0  0  0 -1  1  0  0  0 -1  0  1  0  2 -1  1  0]
        [ 0  0  0  0  0  0 -1  0 -1  0  0  1  0  0 -1  2 -1 -1]
        [ 0  0  0  1  0  0  1  0  0  0  0  0  0 -1  1 -1  2 -1]
        [ 0  0  0 -1  0  0  0  0  1  0  0 -1  0  1  0 -1 -1  2]

    """
    D = incidence_matrix(Gamma, eo)
    return D*D.transpose()
    
def quotient_adjacent(orbits, VG, EG, i, j):
    if (orbits[i],orbits[j]) in EG or (orbits[j],orbits[i]) in EG:
        #print i,j
        return 1
    else:
        return 0


def quotient_graph(Gamma, G, verbose=False):
    """
    This function returns the graph quotient, and in verbose mode
    information on the vertex correspondence of the quotient map.

    If Gamma = (V, E) is a graph and G a subgroup of its 
    automorphism group, we define the quotient graph 
    Gamma/G as follows:
    (1) The vertices of Gamma/G are the G-orbits in V. 
    (2) Distinct vertices v1, v2 of Gamma/G are connected by an 
        edge if and only if there is a vertex v1' in V belonging
        to the orbit of v1, and a vertex v2' in V belonging
        to the orbit of v2, for which (v1', v2') belongs to E.
    (3) Gamma/G is simple.
    
    Input: A graph and a subgroup of its automorphism group.
    Output: The quotient graph.

    EXAMPLES:
        sage: Gamma = graphs.PaleyGraph(13)
        sage: G = Gamma.automorphism_group(); G
        Permutation Group with generators [(1,9,3)(2,5,6)(4,10,12)(7,11,8), (1,10,9,12,3,4)(2,7,5,11,6,8), (0,1,2,3,4,5,6,7,8,9,10,11,12)]
        sage: quotient_graph(Gamma, G)
        Multi-graph on 1 vertex
        sage: G2 = G.sylow_subgroup(2)
        sage: G2.order()
        2
        sage: G3 = G.sylow_subgroup(3)
        sage: G3.order()
        3
        sage: quotient_graph(Gamma, G3)
        Graph on 5 vertices
        sage: quotient_graph(Gamma, G3).show()
        sage: quotient_graph(Gamma, G3, verbose=True)
        (Graph on 5 vertices,
         [0, 1, 2, 3, 4],
         [({8, 11, 7}, {2, 5, 6}),
          ({8, 11, 7}, {1, 3, 9}),
          ({2, 5, 6}, {1, 3, 9}),
          ({0}, {1, 3, 9}),
          ({4, 10, 12}, {8, 11, 7}),
          ({4, 10, 12}, {1, 3, 9}),
          ({4, 10, 12}, {2, 5, 6}),
          ({4, 10, 12}, {0})],
         [{4, 10, 12}, {8, 11, 7}, {2, 5, 6}, {0}, {1, 3, 9}])

    """
    V = Gamma.vertices()
    E = Gamma.edges(labels=False)
    orbits = list(Set([Set([g(v) for g in G]) for v in V]))
    if len(orbits) == 1:
        return Graph(matrix([[0]]))
    n = len(orbits)
    VG = range(n)
    EG = []
    for v in orbits:
        for w in orbits:
            if v<>w:
                for x in v:
                    for y in w:
                        if (x,y) in E or (y,x) in E:
                            if not((w,v) in EG):
                                EG = EG+[(v,w)]
                                break
    EG = list(Set(EG))
    orbs = [list(x) for x in orbits]
    #orbs.sort()
    A = [[quotient_adjacent(orbits, VG, EG, i, j) for i in VG] for j in VG]
    if verbose==True:
        return Graph(matrix(A), format = "adjacency_matrix"), VG, EG, orbs
    else:
        return Graph(matrix(A), format = "adjacency_matrix")


def zeta_function_graph(Gamma):
    """
    Returns the zeta function of the graph Gamma. It is very slow
    except for small graphs.


    EXAMPLES:
        sage: Gamma = graphs.DorogovtsevGoltsevMendesGraph(2)
        sage: zeta_function_graph(Gamma)
        (9*x^4 + 6*x^2 + x + 1)^2*(12*x^3 + 2*x - 1)*(x^2 - 1)^3

    """
    x = var("x")
    A = Gamma.adjacency_matrix()
    r = 1 + len(Gamma.edges()) - len(Gamma.vertices())
    Q = Gamma.laplacian_matrix()
    n = len(Gamma.vertices())
    In = identity_matrix(n)
    return (1-x^2)^(r-1)*(In - A*x + Q*x^2).determinant().factor()

def star_subgraph(Gamma, v):
    """
    Returns the star subgraph of Gamma based at the vertex v.

    EXAMPLES:
        sage: Gamma = graphs.WheelGraph(4)
        sage: Gamma1 = star_subgraph(Gamma, 0)
        sage: Gamma1.vertices()
        [0, 1, 2, 3]
        sage: Gamma1.edges()
        [(0, 1, None), (0, 2, None), (0, 3, None)]

    """
    V = Gamma.vertices()
    E = Gamma.edges()
    Vv = Gamma.neighbors(v)
    Ev = [e for e in E if v in e]
    Ev = [(e[0],e[1]) for e in Ev]
    #print Ev
    return Gamma.subgraph(edges = Ev)

def is_graph_morphism(Gamma1, Gamma2, f):
    """
    Returns True if f defines a graph-theoretic mapping
    from Gamma2 to Gamma 1(ie, vertices to vertices, 
    edges to edges), and False otherwise. 

    Suppose Gamma2 has n vertices and m edges. A morphism 
    f: Gamma2 -> Gamma1
    can represented by 2 pairs of lists [VL2, VL1] and [EL2, EL1],
    where VL2 is the list of all n vertices of Gamma2,
    VL1 is the list of length n of the vertices
    in Gamma1 that form the corresponding image under
    the map f, EL2 is the list of all m edges of Gamma2,
    EL1 is the list of length m of the edges
    in Gamma1 that form the corresponding image under
    the map f (including None). However, the pair [EL2, EL1]
    can be deduced from the pair [VL2, VL1], so we 
    drop [EL2, EL1].

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: is_graph_morphism(Gamma1, Gamma2, f)
        True
        sage: Gamma2 = graphs.WheelGraph(5)
        sage: Gamma1 = graphs.WheelGraph(2)
        sage: f = [[0,1,2,3,4],[0,1,1,1,1]]
        sage: is_graph_morphism(Gamma1, Gamma2, f)
        True

    """
    V1 = Gamma1.vertices()
    V2 = Gamma2.vertices()
    E1 = Gamma1.edges()
    EE1 = [(x[0],x[1]) for x in E1]
    E2 = Gamma2.edges()
    domainfV = f[0]
    imagefV = f[1]
    for i in domainfV:
        if not(i in V2):
            print "Not a function on all of Gamma2"
            return False
    for i in imagefV:
        if not(i in V1):
            print "Not a function whose image is in Gamma1"
            return False
    for e in E2:
        x2 = e[0]
        y2 = e[1]
        i2 = V2.index(x2)
        j2 = V2.index(y2)
        i1 = imagefV[i2]
        j1 = imagefV[j2]
        if i1<>j1 and not((i1,j1) in EE1) and not((j1,i1) in EE1):
            print (i1,j1),  EE1
            print "Does not send edges to {edges} \cup {vertices}"
            return False
    return True

def edge_image_list(Gamma1, Gamma2, f):
    """

    EXAMPLES:
        sage: Gamma2 = graphs.WheelGraph(5)
        sage: Gamma1 = graphs.WheelGraph(2)
        sage: f = [[0,1,2,3,4],[0,1,1,1,1]]
        sage: is_graph_morphism(Gamma1, Gamma2, f)
        True
        sage: edge_image_list(Gamma1, Gamma2, f)
        [[(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (3, 4)],
         [(0, 1), (0, 1), (0, 1), (0, 1), 'None', 'None', 'None', 'None']]

    """
    V1 = Gamma1.vertices()
    V2 = Gamma2.vertices()
    E1 = Gamma1.edges()
    EE1 = [(x[0],x[1]) for x in E1]
    E2 = Gamma2.edges()
    EE2 = [(x[0],x[1]) for x in E2]
    domainfV = f[0]
    imagefV = f[1]
    imageE2 = []
    for e in E2:
        x2 = e[0]
        y2 = e[1]
        i2 = V2.index(x2)
        j2 = V2.index(y2)
        i1 = imagefV[i2]
        j1 = imagefV[j2]
        if (i1,j1) in EE1 or (j1,i1) in EE1:
            imageE2.append((i1,j1))
        else:
            imageE2.append("None")
    return [EE2, imageE2]

#def image_of_edges(Gamma1, Gamma2, f):
#    L = []
#    for e in Gamma2.edges():
#        L.append(image_of_edge_under_graph_morphism(Gamma1, Gamma2, f, e))
#    return L

def index_of_graph_morphism_edges(Gamma1, Gamma2, f, e2, e1, e0):
    """
    These are the entries of the signed edge matrix attached to the 
    harmonic map f : Gamma2 -> Gamma1.

    Here e1 is an edge in Gamma1 and e2 is an edge in Gamma2,
    and e0 is an orientation on Gamma1.
    """
    Ef = edge_image_list(Gamma1, Gamma2, f)
    V1 = Gamma1.vertices()
    E1 = Gamma1.edges()
    EE1 = [(x[0],x[1]) for x in E1]
    i_e1 = V1.index(e1[0])
    j_e1 = EE1.index(e1)
    if e2 in Ef[0]:
        i = Ef[0].index(e2)
    else:
        i = Ef[0].index((e2[1],e2[0]))
    if not(e2 in Ef[0]):
        raise ValueError("argument is not an edge of the graph")
    if (e1[0],e1[1]) == Ef[1][i]:
        return e0[i_e1]
    elif (e1[1],e1[0]) == Ef[1][i]:
        return -e0[i_e1]
    else:
        return 0

def index_of_graph_morphism_vertices(Gamma1, Gamma2, f, v2, v1):
    i =f[0].index(v2)
    if v1 == f[1][i]:
        return 1
    else:
        return 0

def matrix_of_graph_morphism_vertices(Gamma1, Gamma2, f):
    """
    Returns an n2 x n1 (0,1)-matrix representing Phi_V 
    if f defines a graph-theoretic mapping
    from Gamma2 to Gamma 1 (ie, vertices to vertices, 
    edges to edges), and False otherwise. 

    Suppose Gamma2 has n vertices and m edges. A morphism 
    f: Gamma2 -> Gamma1
    can represented by 2 pairs of lists [VL2, VL1] and [EL2, EL1],
    where VL2 is the list of all n vertices of Gamma2,
    VL1 is the list of length n of the vertices
    in Gamma1 that form the corresponding image under
    the map f, EL2 is the list of all m edges of Gamma2,
    EL1 is the list of length m of the edges
    in Gamma1 that form the corresponding image under
    the map f (including None). However, the pair [EL2, EL1]
    can be deduced from the pair [VL2, VL1].

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: is_graph_morphism(Gamma1, Gamma2, f)
        True
        sage: Phi = matrix_of_graph_morphism_vertices(Gamma1, Gamma2, f); Phi
        [1 0 1 0]
        [0 1 0 1]
        sage: Phi.transpose()*Phi
        [1 0 1 0]
        [0 1 0 1]
        [1 0 1 0]
        [0 1 0 1]
        sage: Phi*Phi.transpose()
        [2 0]
        [0 2]
        sage: Gamma2 = graphs.CubeGraph(3)
        sage: Gamma1 = graphs.TetrahedralGraph()
        sage: f = [['000', '001', '010', '011', '100', '101', '110', '111'], [0, 1, 2, 3, 3, 2, 1, 0]]
        sage: Phi = matrix_of_graph_morphism_vertices(Gamma1, Gamma2, f); Phi
        [1 0 0 0 0 0 0 1]
        [0 1 0 0 0 0 1 0]
        [0 0 1 0 0 1 0 0]
        [0 0 0 1 1 0 0 0]
        sage: D = vector([-1,0,1,2,2,1,0,-1])
        sage: Phi*D
        (-2, 0, 2, 4)
        sage: Phi*Phi.transpose()
        [2 0 0 0]
        [0 2 0 0]
        [0 0 2 0]
        [0 0 0 2]
        sage: Phi.transpose()*Phi
        [1 0 0 0 0 0 0 1]
        [0 1 0 0 0 0 1 0]
        [0 0 1 0 0 1 0 0]
        [0 0 0 1 1 0 0 0]
        [0 0 0 1 1 0 0 0]
        [0 0 1 0 0 1 0 0]
        [0 1 0 0 0 0 1 0]
        [1 0 0 0 0 0 0 1]
        sage: Gamma2 = graphs.WheelGraph(5)
        sage: Gamma1 = graphs.WheelGraph(2)
        sage: f = [[0,1,2,3,4],[0,1,1,1,1]]
        sage: Phi = matrix_of_graph_morphism_vertices(Gamma1, Gamma2, f); Phi
        [1 0 0 0 0]
        [0 1 1 1 1]
        sage: Phi.transpose()*Phi
        [1 0 0 0 0]
        [0 1 1 1 1]
        [0 1 1 1 1]
        [0 1 1 1 1]
        [0 1 1 1 1]
        sage: Phi*Phi.transpose()
        [1 0]
        [0 4]
        sage: Q1 = Gamma1.laplacian_matrix()
        sage: Q2 = Gamma2.laplacian_matrix()
        sage: Q1*Phi
        [ 1 -1 -1 -1 -1]
        [-1  1  1  1  1]
        sage: Phi*Q2
        [ 4 -1 -1 -1 -1]
        [-4  1  1  1  1]
        sage: Gamma3 = graphs.CycleGraph(4)
        sage: Gamma3.add_edge([0,2])
        sage: f2 = [[0,1,2,3,4],[0,2,1,2,3]]
        sage: Phi2 = matrix_of_graph_morphism_vertices(Gamma3, Gamma2, f2); Phi2
        [1 0 0 0 0]
        [0 0 1 0 0]
        [0 1 0 1 0]
        [0 0 0 0 1]
        sage: Phi2.transpose()*Phi2
        [1 0 0 0 0]
        [0 1 0 1 0]
        [0 0 1 0 0]
        [0 1 0 1 0]
        [0 0 0 0 1]
        sage: Phi2*Phi2.transpose()
        [1 0 0 0]
        [0 1 0 0]
        [0 0 2 0]
        [0 0 0 1]
        sage: v = Q3.right_eigenspaces()[0][1].basis()[0]; v
        (0, 1, 0, -1)
        sage: Q3*v ## Does Phi2 send eigenspaces to eigenspaces?
        (0, 2, 0, -2)
        sage: Phi2.transpose()*Phi2*Q2*Phi2.transpose()*v
        (0, 0, 3, 0, -3)
        sage: w = Q3.right_eigenspaces()[2][1].basis()[0]; w
        (1, 0, -1, 0)
        sage: Q3*w # A: No.
        (4, 0, -4, 0)
        sage: Phi2.transpose()*Phi2*Q2*Phi2.transpose()*w
        (6, -8, 1, -8, 1)

    """
    if not(is_graph_morphism(Gamma1, Gamma2, f)):
         raise ValueError("the map given is not a graph morphism")
    V1 = Gamma1.vertices()
    V2 = Gamma2.vertices()
    n1 = len(V1)
    n2 = len(V2)
    MS = MatrixSpace(ZZ, n1, n2)
    L = [[index_of_graph_morphism_vertices(Gamma1, Gamma2, f, v2, v1) for v2 in V2] for v1 in V1]
    return MS(L).transpose()


def matrix_of_graph_morphism_edges(Gamma1, Gamma2, f, e0):
    """
    Returns an m2 x m1 (0,1)-matrix representing Phi_E 
    if f defines a graph-theoretic mapping
    from Gamma2 to Gamma 1 (ie, vertices to vertices, 
    edges to edges), and False otherwise. 

    Here e0 is an orientation on Gamma1. 

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: eo = [1]
        sage: is_graph_morphism(Gamma1, Gamma2, f)
        True
        sage: PhiE = matrix_of_graph_morphism_edges(Gamma1, Gamma2, f, eo); PhiE
        [1]
        [0]
        [0]
        [1]
        sage: PhiE.transpose()*PhiE
        [2]

        sage: Gamma2 = graphs.CubeGraph(3)
        sage: Gamma1 = graphs.TetrahedralGraph()
        sage: eo = [1,1,1,1,1,1]
        sage: f = [['000', '001', '010', '011', '100', '101', '110', '111'], [0, 1, 2, 3, 3, 2, 1, 0]]
        sage: PhiE = matrix_of_graph_morphism_edges(Gamma1, Gamma2, f, eo); PhiE
        [ 1  0  0  0  0  0]
        [ 0  1  0  0  0  0]
        [ 0  0  1  0  0  0]
        [ 0  0  0  0  1  0]
        [ 0  0  0  1  0  0]
        [ 0  0  0  0  0  1]
        [ 0  0  0 -1  0  0]
        [ 0  0 -1  0  0  0]
        [ 0  0  0  0  0 -1]
        [ 0  0  0  0 -1  0]
        [ 0 -1  0  0  0  0]
        [-1  0  0  0  0  0]
        sage: PhiE.transpose()*PhiE
        [2 0 0 0 0 0]
        [0 2 0 0 0 0]
        [0 0 2 0 0 0]
        [0 0 0 2 0 0]
        [0 0 0 0 2 0]
        [0 0 0 0 0 2]

        sage: Gamma2 = graphs.WheelGraph(5)
        sage: Gamma1 = graphs.WheelGraph(2)
        sage: f = [[0,1,2,3,4],[0,1,1,1,1]]
        sage: eo = [1]
        sage: PhiE = matrix_of_graph_morphism_edges(Gamma1, Gamma2, f, eo); PhiE.transpose()
        [1 1 1 1 0 0 0 0]
 
        sage: Gamma3 = graphs.CycleGraph(4)
        sage: Gamma2 = graphs.WheelGraph(5)
        sage: Gamma3.add_edge([0,2])
        sage: f2 = [[0,1,2,3,4],[0,2,1,2,3]]
        sage: eo = [1]*5
        sage: PhiE2 = matrix_of_graph_morphism_edges(Gamma3, Gamma2, f2, eo); PhiE2.transpose()
        [ 0 -1  0  0  0  0  0  0]
        [-1  0 -1  0  0  0  0  0]
        [ 0  0  0 -1  0  0  0  0]
        [ 0  0  0  0 -1  0 -1  0]
        [ 0  0  0  0  0 -1  0 -1]

        sage: Gamma2 = graphs.PaleyGraph(9); Gamma2
        Paley graph with parameter 9: Graph on 9 vertices
        sage: Gamma1 = graphs.CycleGraph(3); Gamma1
        Cycle graph: Graph on 3 vertices
        sage: eo = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        sage: eoo = [1,1,1]
        sage: B1 = incidence_matrix(Gamma1, eoo)
        sage: B2 = incidence_matrix(Gamma2, eo)
        sage: V2 = Gamma2.vertices()
        sage: f = [V2, [0,0,0,2,2,2,1,1,1]]
        sage: is_harmonic_graph_morphism(Gamma1, Gamma2, f) 
        True
        sage: PhiV = matrix_of_graph_morphism_vertices(Gamma1, Gamma2, f)
        sage: PhiE = matrix_of_graph_morphism_edges(Gamma1, Gamma2, f, eoo)
        sage: B1*PhiV.transpose()
        [ 1  1  1  0  0  0 -1 -1 -1]
        [ 1  1  1 -1 -1 -1  0  0  0]
        [ 0  0  0 -1 -1 -1  1  1  1]
        sage: PhiE.transpose()*B2
        [ 1  1  1  0  0  0 -1 -1 -1]
        [ 1  1  1 -1 -1 -1  0  0  0]
        [ 0  0  0 -1 -1 -1  1  1  1]

    """
    if not(is_graph_morphism(Gamma1, Gamma2, f)):
         raise ValueError("the map given is not a graph morphism")
    E1 = Gamma1.edges()
    E2 = Gamma2.edges()
    EE1 = [(x[0],x[1]) for x in E1]
    EE2 = [(x[0],x[1]) for x in E2]
    m1 = len(E1)
    m2 = len(E2)
    MS = MatrixSpace(ZZ, m1, m2)
    L = [[index_of_graph_morphism_edges(Gamma1, Gamma2, f, e2, e1, e0) for e2 in EE2] for e1 in EE1]
    return MS(L).transpose()

def preimage_of_vertex_under_graph_morphism(Gamma1, Gamma2, f, v):
    """
    Returns the preimage(s) of v under f:Gamma2 -> Gamma 1,
    where v is an vertex in Gamma1.

    A morphism  f: Gamma2 -> Gamma1
    is represented by a pair of lists [L2, L1],
    where L2 is the list of all n vertices of Gamma2,
    and L1 is the list of length n of the vertices
    in Gamma1 that form the corresponding image under
    the map f.

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: Gamma1.vertices()
        [('00', '01', None)]
        sage: v = Gamma1.vertices()[0]
        sage: preimage_of_vertex_under_graph_morphism(Gamma1, Gamma2, f, v)
        ['00', '10']

    """
    V1 = Gamma1.vertices()
    V2 = Gamma2.vertices()
    domainf = f[0]
    imagef = f[1]
    indices = [i for i in range(len(V2)) if imagef[i]==v]
    return [V2[i] for i in indices]

def preimage_of_edge_under_graph_morphism(Gamma1, Gamma2, f, e):
    """
    Returns the preimage(s) of e under f:Gamma2 -> Gamma 1,
    where e is an edge in Gamma1.

    A morphism  f: Gamma2 -> Gamma1
    is represented by a pair of lists [L2, L1],
    where L2 is the list of all n vertices of Gamma2,
    and L1 is the list of length n of the vertices
    in Gamma1 that form the corresponding image under
    the map f.

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1.edges()
        [('00', '01', None)]
        sage: e = Gamma1.edges()[0]
        sage: preimage_of_edge_under_graph_morphism(Gamma1, Gamma2, f, e)
        [('00', '01', None), ('10', '11', None)]
        sage: Gamma2 = graphs.CubeGraph(3)
        sage: Gamma1 = graphs.TetrahedralGraph()
        sage: f = [['000', '001', '010', '011', '100', '101', '110', '111'], [0, 1, 2, 3, 3, 2, 1, 0]]
        sage: e = Gamma1.edges()[1]; e
        (0, 2, None)
        sage: preimage_of_edge_under_graph_morphism(Gamma1, Gamma2, f, e)
        [('000', '010', None), ('101', '111', None)]

    """
    V1 = Gamma1.vertices()
    V2 = Gamma2.vertices()
    E1 = Gamma1.edges()
    E2 = Gamma2.edges()
    domainf = f[0]
    imagef = f[1]
    imagef_e = []
    vs = preimage_of_vertex_under_graph_morphism(Gamma1, Gamma2, f, e[0])
    ws = preimage_of_vertex_under_graph_morphism(Gamma1, Gamma2, f, e[1])
    edges = []
    for x in vs:
        for y in ws:
            for e in E2:
                if (x,y) == (e[0],e[1]) or (x,y) == (e[1],e[0]):
                    edges.append(e)
    return edges

def image_of_edge_under_graph_morphism(Gamma1, Gamma2, f, e):
    """
    Returns the image of e under f:Gamma2 -> Gamma 1,
    where e is an edge in Gamma2.

    A morphism  f: Gamma2 -> Gamma1
    is represented by a pair of lists [L2, L1],
    where L2 is the list of all n vertices of Gamma2,
    and L1 is the list of length n of the vertices
    in Gamma1 that form the corresponding image under
    the map f.

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1.edges()
        [('00', '01', None)]
        sage: e = Gamma1.edges()[0]; e
        ('00', '01', None)
        sage: image_of_edge_under_graph_morphism(Gamma1, Gamma2, f, e)
        ('00', '01', None)
        sage: e = Gamma2.edges()[1]; e
        ('00', '10', None)
        sage: image_of_edge_under_graph_morphism(Gamma1, Gamma2, f, e)
        ('00', '00')


    """
    V1 = Gamma1.vertices()
    V2 = Gamma2.vertices()
    E1 = Gamma1.edges()
    E2 = Gamma2.edges()
    domainf = f[0]
    imagef = f[1]
    imagef_e = []
    x2 = e[0]
    y2 = e[1]
    i2 = V2.index(x2)
    j2 = V2.index(y2)
    i1 = imagef[i2]
    j1 = imagef[j2]
    for e2 in E1:
        if (i1,j1) == (e2[0],e2[1]) or (j1,i1) == (e2[1],e2[0]):
            return e
    return (i1,j1)

def image_of_vertex_under_graph_morphism(Gamma1, Gamma2, f, v):
    """
    Returns the image of v under f:Gamma2 -> Gamma 1,
    where v is an vertex in Gamma2.

    A morphism  f: Gamma2 -> Gamma1
    is represented by a pair of lists [L2, L1],
    where L2 is the list of all n vertices of Gamma2,
    and L1 is the list of length n of the vertices
    in Gamma1 that form the corresponding image under
    the map f.

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: v = Gamma2.vertices()[0]
        sage: image_of_vertex_under_graph_morphism(Gamma1, Gamma2, f, v)
        '00'

    """
    V1 = Gamma1.vertices()
    V2 = Gamma2.vertices()
    E1 = Gamma1.edges()
    E2 = Gamma2.edges()
    domainf = f[0]
    imagef = f[1]
    i = V2.index(v)
    j = imagef[i]
    return j

def is_harmonic_graph_morphism(Gamma1, Gamma2, f, verbose = False):
    """
    Returns True if f defines a graph-theoretic mapping
    from Gamma2 to Gamma1 that is harmonic, and False otherwise. 

    Suppose Gamma2 has n vertices. A morphism 
              f: Gamma2 -> Gamma1
    is represented by a pair of lists [L2, L1],
    where L2 is the list of all n vertices of Gamma2,
    and L1 is the list of length n of the vertices
    in Gamma1 that form the corresponding image under
    the map f.

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: is_harmonic_graph_morphism(Gamma1, Gamma2, f)
        True
        sage: Gamma2 = graphs.CubeGraph(3)
        sage: Gamma1 = graphs.TetrahedralGraph()
        sage: f = [['000', '001', '010', '011', '100', '101', '110', '111'], [0, 1, 2, 3, 3, 2, 1, 0]]
        sage: is_harmonic_graph_morphism(Gamma1, Gamma2, f)
        True
        sage: Gamma2 = graphs.CubeGraph(3)
        sage: Gamma1 = graphs.CubeGraph(2)
        sage: f = [['000', '001', '010', '011', '100', '101', '110', '111'], ["00", "00", "01", "01", "10", "10", "11", "11"]]
        sage: is_harmonic_graph_morphism(Gamma1, Gamma2, f)
        True
        sage: is_harmonic_graph_morphism(Gamma1, Gamma2, f, verbose=True)
        This [<vertex>, <intersection sizes>]] passes the check: ['000', [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: ['001', [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: ['010', [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: ['011', [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: ['100', [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: ['101', [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: ['110', [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: ['111', [1, 1]]
        True
        sage: Gamma2 = graphs.TetrahedralGraph()
        sage: Gamma1 = graphs.CycleGraph(3)
        sage: f = [[0,1,2,3],[0,1,2,0]]
        sage: is_harmonic_graph_morphism(Gamma1, Gamma2, f)
        False
        sage: is_harmonic_graph_morphism(Gamma1, Gamma2, f, verbose=True)
        This [<vertex>, <intersection sizes>]] passes the check: [0, [1, 1]]
        This [<vertex>, <intersection sizes>]] fails the check: [1, [2, 1]]
        This [<vertex>, <intersection sizes>]] fails the check: [2, [2, 1]]
        False
        sage: V1 = [0,1,2]
        sage: V2 = [0,1,2,3]
        sage: E2 = [(0,1),(1,2),(2,0),(2,3)]
        sage: E1 = [(0,1),(1,2),(2,0)]
        sage: Gamma2 = Graph([V2,E2],format='vertices_and_edges')
        sage: Gamma1 = Graph([V1,E1],format='vertices_and_edges')
        sage: f = [V2, [0,1,2,2]]    # thi is a retract
        sage: is_harmonic_graph_morphism(Gamma1, Gamma2, f, verbose=True)
        This [<vertex>, <intersection sizes>]] passes the check: [0, [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: [1, [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: [2, [1, 1]]
        This [<vertex>, <intersection sizes>]] passes the check: [3, [0, 0]]
        True

    """
    V1 = Gamma1.vertices()
    n1 = len(V1)
    V2 = Gamma2.vertices()
    n2 = len(V2)
    E1 = Gamma1.edges()
    m1 = len(E1)
    E2 = Gamma2.edges()
    m2 = len(E2)
    edges_in_common = []
    for v2 in V2:
        w = image_of_vertex_under_graph_morphism(Gamma1, Gamma2, f, v2)
        str1 = star_subgraph(Gamma1, w)
        Ew = str1.edges()
        str2 = star_subgraph(Gamma2, v2)
        Ev2 = str2.edges()
        sizes = []
        for e in Ew:
            finv_e = preimage_of_edge_under_graph_morphism(Gamma1, Gamma2, f, e)
            L = [x for x in finv_e if x in Ev2]
            sizes.append(len(L))
            #print v2,e,L
        edges_in_common.append([v2, sizes]) # these common values = horiz. mult.
    ans = True
    for x in edges_in_common:
        sizes = x[1]
        S = Set(sizes)
        if S.cardinality()>1:
            ans = False
            if verbose and ans==False:
                print "This [<vertex>, <intersection sizes>]] fails the check:", x
        if verbose and ans==True:
            print "This [<vertex>, <intersection sizes>]] passes the check:", x
    return ans
            
def horizonal_multiplicity_of_harmonic_map(Gamma1, Gamma2, f, v2):
    """
    Returns m_\phi(v)
    """
    V1 = Gamma1.vertices()
    n1 = len(V1)
    V2 = Gamma2.vertices()
    n2 = len(V2)
    E1 = Gamma1.edges()
    m1 = len(E1)
    E2 = Gamma2.edges()
    m2 = len(E2)
    edges_in_common = []
    w = image_of_vertex_under_graph_morphism(Gamma1, Gamma2, f, v2)
    str1 = star_subgraph(Gamma1, w)
    Ew = str1.edges()
    str2 = star_subgraph(Gamma2, v2)
    Ev2 = str2.edges()
    sizes = []
    for e in Ew:
        finv_e = preimage_of_edge_under_graph_morphism(Gamma1, Gamma2, f, e)
        L = [x for x in finv_e if x in Ev2]
        sizes.append(len(L))
        #print v2,e,L
    edges_in_common.append([v2, sizes]) # these common values = horiz. mult.
    return edges_in_common[0][1][0]

def genus(Gamma):
    """
    Not the same as Gamma.genus()

    EXAMPLES:
        sage: Gamma2 = graphs.WheelGraph(5)
        sage: Gamma2.genus()
        0
        sage: genus(Gamma2)
        4

    """
    return len(Gamma.edges()) - len(Gamma.vertices()) + 1

def perm_number_action(G, g0):
    Glist = G.list()
    N = len(Glist)
    L = [(Glist.index(x),Glist.index(x*g0)) for x in G]
    return L

def lifted_vertex_to_cover(G, Gamma, v, g0):
    """
    In a (to be constructed) G-covering of Gamma, 
    call it Gamma' = (V', E'), this function 
    returns the vertex in V' that v lifts to associated to
    g0 under the inverse of the projection map 
    phi: Gamma' -> Gamma.

    INPUT:
        G - a non-trivial subgroup of a symmetric group,
        Gamma - a graph with vertices V and edges E
        v - a vertex in V
        g0 - a non-trivial element of G

    EXAMPLES:

    """
    V1 = Gamma.vertices()
    Glist = G.list()
    n = G.cardinality()
    V2 = range(n*len(V1))
    iv = V1.index(v)
    j = perm_number_action(G, g0)[iv][1]
    return iv+j


def lifted_vertices_to_cover(G, Gamma, v, g0 = 1):
    """
    In a (to be constructed) G-covering of Gamma, 
    call it Gamma' = (V', E'), this function 
    returns the vertex in V' that v lifts to associated to
    the inverse of the projection map 
    phi: Gamma' -> Gamma.

    INPUT:
        G - a non-trivial subgroup of a symmetric group,
        Gamma - a graph with vertices V and edges E
        v - a vertex in V
        g0 - a non-trivial element of G

    EXAMPLES:

    """
    V1 = Gamma.vertices()
    Glist = G.list()
    n = G.cardinality()
    n1 = len(V1)
    V2 = range(n*n1)
    iv = V1.index(v)
    if g0<>1:
        iv = g0(iv)
        #print iv
    return [(iv+j*n1)%(n*n1) for j in range(n)]

def lifted_edges_to_cover(G, Gamma, e):
    """
    In a (to be constructed) G-covering of Gamma, 
    call it Gamma' = (V', E'), this function 
    returns the edges in E' that e lifts to associated to
    the inverse of the projection map 
    phi: Gamma' -> Gamma.

    INPUT:
        G - a non-trivial subgroup of a symmetric group,
        Gamma - a graph with vertices V and edges E
        e - a vertex in E

    EXAMPLES:

    """
    v = e[0]
    w = e[1]
    V1 = Gamma.vertices()
    Glist = G.list()
    n = G.cardinality()
    n1 = len(V1)
    cyclics_Sn1 = [g0 for g0 in SymmetricGroup(n1) if g0.order()==3]
    Lv = lifted_vertices_to_cover(G, Gamma, v)
    Lw = lifted_vertices_to_cover(G, Gamma, w, cyclics_Sn1[0])
    #print v,w,Lv, Lw, g00
    E2 = []
    for i in range(n):
        if Lv[i]<>Lw[i] and not((Lv[i],Lw[i]) in E2) and not((Lw[i],Lv[i]) in E2):
            E2.append([Lv[i],Lw[i]])
    Lv = lifted_vertices_to_cover(G, Gamma, v, cyclics_Sn1[0]^(-1))
    Lw = lifted_vertices_to_cover(G, Gamma, w)
    for i in range(n):
        if Lv[i]<>Lw[i] and not((Lv[i],Lw[i]) in E2) and not((Lw[i],Lv[i]) in E2):
            E2.append([Lv[i],Lw[i]])
    return E2

def lifted_graph_to_cover0(G, Gamma):
    """
    Constructs a G-covering of Gamma.

    INPUT:
        G - a non-trivial subgroup of a symmetric group,
        Gamma - a graph with vertices V and edges E

    EXAMPLES:
        sage: Gamma = graphs.CycleGraph(3)
        sage: S3 = SymmetricGroup(3)
        sage: Gamma2 = lifted_graph_to_cover0(S3, Gamma)
        sage: Gamma2.show()
        sage: Gamma2.automorphism_group().cardinality()
        768
        sage: quotient_graph(Gamma2, AG.sylow_subgroup(3))
        Graph on 6 vertices

    """
    V1 = Gamma.vertices()
    n = G.cardinality()
    n1 = len(V1)
    Gamma2 = Graph(n*n1)
    E2 = []
    for e in Gamma.edges():
        e = [e[0],e[1]]
        es = lifted_edges_to_cover(G, Gamma, e)
        #print e, len(es)
        E2 = E2+es
    Gamma2.add_edges(E2)
    #print E2, Gamma2.vertices()
    return Gamma2


def lifted_graph_to_cover(G, Gamma):
    """
    Constructs a G-covering of Gamma.

    INPUT:
        G - a non-trivial subgroup of a symmetric group,
        Gamma - a graph with vertices V and edges E

    EXAMPLES:
        sage: Gamma = graphs.CycleGraph(3)
        sage: S3 = SymmetricGroup(3)
        sage: Gamma2 = lifted_graph_to_cover(S3, Gamma)
        sage: Gamma2.show()
        sage: Gamma2.automorphism_group().cardinality()
        768
        sage: quotient_graph(Gamma2, AG.sylow_subgroup(3))
        Graph on 6 vertices

    """
    V1 = Gamma.vertices()
    n = G.cardinality()
    n1 = len(V1)
    Gamma2 = Graph(n*n1)
    E2 = []
    for e in Gamma.edges():
        e = [e[0],e[1]]
        es = lifted_edges_to_cover(G, Gamma, e)
        #print e, len(es)
        E2 = E2+es
    Gamma2.add_edges(E2)
    #print E2, Gamma2.vertices()
    return Gamma2



###################################################
####### Cayley graph for linear codes
###################################################

def cayley_expander_graph_from_gen_mat(G):
    """

    EXAMPLES:
        sage: C = codes.HammingCode(3, GF(2)); C
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: G = C.generator_matrix()
        sage: cayley_expander_graph_from_gen_mat(G)
        Graph on 16 vertices

    """
    Gcols = G.columns()
    n = len(Gcols)
    k = len(G.rows()) # assumes G is full rank
    V = GF(2)**k
    f = lambda x: ZZ(x in Gcols)
    C = cartesian_product([range(2**k), range(2**k)])
    Vlist = V.list()                  
    E = [(x[0],x[1]) for x in C if f(Vlist[x[0]]+Vlist[x[1]])==1]
    E = Set([Set(s) for s in E])
    E = [tuple(s) for s in E] 
    Gamma = Graph(E)
    return Gamma

def cayley_expander_graph_from_code(C):
    """

    EXAMPLES:
       sage: C = codes.HammingCode(3, GF(2)); C
       Linear code of length 7, dimension 4 over Finite Field of size 2
       sage: Gamma = cayley_expander_graph_from_code(C)
       sage: Gamma.show(layout="circular")
       Launched png viewer for Graphics object consisting of 73 graphics primitives


    """
    G = C.generator_matrix()
    Gcols = G.columns()
    n = C.length()
    k = C.dimension()
    f = lambda x: ZZ(x in G.columns())
    V = GF(2)**k
    CP = cartesian_product([range(2**k), range(2**k)])
    Vlist = V.list()                  
    E = [(x[0],x[1]) for x in CP if f(Vlist[x[0]]+Vlist[x[1]])==1]
    E = Set([Set(s) for s in E])
    E = [tuple(s) for s in E] 
    Gamma = Graph(E)
    return Gamma

def edge_expansion(Gamma):
    """
    """
    V = Gamma.vertices()
    SV = Set(V) 
    return min([(QQ(len(Gamma.edge_boundary([V[i] for i in SV.list()])))/QQ(len(SV.list()))) for SV in SetV.subsets() if (len(SV.list())<>0 and len(SV.list())<5)])

def laplacian_pseudoinverse(Gamma):
    """
    Returns Q^+

    EXAMPLES:
        sage: Gamma2 = graphs.CycleGraph(4); Gamma2
        Cycle graph: Graph on 4 vertices
        sage: latex(Gamma2.laplacian_matrix())
        \left(\begin{array}{rrrr}
        2 & -1 & 0 & -1 \\
        -1 & 2 & -1 & 0 \\
        0 & -1 & 2 & -1 \\
        -1 & 0 & -1 & 2
        \end{array}\right)
        sage: latex(laplacian_pseudoinverse(Gamma2))
        \left(\begin{array}{rrrr}
        \frac{5}{16} & -\frac{1}{16} & -\frac{3}{16} & -\frac{1}{16} \\
        -\frac{1}{16} & \frac{5}{16} & -\frac{1}{16} & -\frac{3}{16} \\
        -\frac{3}{16} & -\frac{1}{16} & \frac{5}{16} & -\frac{1}{16} \\
        -\frac{1}{16} & -\frac{3}{16} & -\frac{1}{16} & \frac{5}{16}
        \end{array}\right)

    """
    Q = Gamma.laplacian_matrix()
    n = len(Q.rows())
    J = matrix(QQ, n, n, n^2*[1])
    Qplus = (Q+(1/n)*J)^(-1)-(1/n)*J
    return Qplus

def hamiltonian_paths0(Gamma):
    """
    Returns a list of hamiltonian paths (spanning trees of 
    max degree <=2).

    EXAMPLES:
        sage: Gamma = graphs.GridGraph([3,3])
        sage: HP = hamiltonian_paths(Gamma)
        sage: len(HP)
        20
    """
    ST = Gamma.spanning_trees()
    #print len(ST)
    HP = []
    for X in ST:
        L = X.degree_sequence()
        #print L,ST.index(X), max(L)
        if max(L)<=2:
            #print L,ST.index(X), max(L)
            HP.append(X)
    return HP
    
def hamiltonian_paths(Gamma, signed_adjacency_matrix = []):
    """
    Returns a list of hamiltonian paths (spanning trees of 
    max degree <=2).

    EXAMPLES:
        sage: Gamma = graphs.GridGraph([3,3])
        sage: HP = hamiltonian_paths(Gamma)
        sage: len(HP)
        20
        sage: A = matrix(QQ,[
        [0 , -1 , 1  , -1 , -1 , -1 ],
        [1,   0 ,  -1,  1,  1,   -1  ],
        [-1 , 1 ,  0 ,  1 , 1  , -1  ],
        [1 , -1 , -1,  0 ,  -1 , -1  ],
        [1 , - 1 , - 1 , 1 , 0 , - 1  ],
        [1 ,  1  ,  1  , 1  , 1  , 0 ]
        ])
        sage: Gamma = Graph(A, format='weighted_adjacency_matrix')
        sage: HP = hamiltonian_paths(Gamma, signed_adjacency_matrix = A)
        sage: L = [sum(x[2]) for x in HP]; max(L)
        5
        sage: L.index(5)
        21
        sage: HP[21]                                 
        [Graph on 6 vertices,
         [0, 5, 2, 1, 3, 4],
         [-1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1]]
        sage: L.count(5)
        1

    """
    ST = Gamma.spanning_trees()
    if signed_adjacency_matrix == []:
        HP = []
        for X in ST:
            L = X.degree_sequence()
            if max(L)<=2:
                #print L,ST.index(X), max(L)
                HP.append(X)
        return HP
    if signed_adjacency_matrix != []:
        A = signed_adjacency_matrix
        HP = []
        for X in ST:
            L = X.degree_sequence()
            if max(L)<=2:
                #VX = X.vertices()
                EX = X.edges()
		if EX[0][1] != EX[-1][1]:
                    ranking = X.shortest_path(EX[0][0],EX[-1][1])
		else:
		    ranking = X.shortest_path(EX[0][0],EX[-1][0])
		signature = [A[ranking[i]][ranking[j]] for i in range(len(ranking)-1) for j in range(i+1,len(ranking))]
		#print X,ranking,signature
		# (use X.shortest_path(a,b) instead of path)
                HP.append([X,ranking,signature])
        return HP

def minimum_upset(M):
    """
    EXAMPLES:
        sage: M = matrix(QQ,[
        [0 , 0 , 1  , 0 , 0 , 0 ],
        [1,   0 ,  0,  1,  1,   0  ],
        [0 , 1 ,  0 ,  1 , 1  , 0  ],
        [1 , 0 , 0,  0 ,  0 , 0  ],
        [1 , 0 , 0 , 1 , 0 , 0  ],
        [1 ,  1  ,  1  , 1  , 1  , 0 ]
        ])
        sage: minimum_upset(M)
        (
        [1 1 0 1 1 1]                    
        [0 1 0 1 0 1]                    
        [0 0 0 0 1 0]                    
        [1 1 0 0 0 1]                    
        [1 0 0 0 0 1]                    
        [1 0 0 0 0 0], [5, 1, 2, 3, 0, 4]
        )
    """
    n = len(M.rows())
    Sn = SymmetricGroup(n)
    # rerank by win-loss
    M1 = M
    def cmp1(L1,L2):
       int(bool(sum(L1)>sum(L2)))
    rowsM = [list(x) for x in M1.rows()]
    rowsM.sort()
    rowsM.reverse()
    M1 = matrix(QQ,rowsM)
    #print M1,"\n"
    wins = sum([sum([M1[j][i] for i in range(j,n)]) for j in range(n)])
    total_wins = sum([sum([M1[j][i] for i in range(n)]) for j in range(n)])
    p = Sn(1)
    for g in Sn:
        P = g.matrix()
        M0 = P*M1*P^(-1)
        if sum([sum([M0[j][i] for i in range(j,n)]) for j in range(n)])>wins:
            if sum([sum([M0[j][i] for i in range(n)]) for j in range(n)])>=total_wins:
                M1 = M0
                p = p*g
    p0 = g^(-1)
    return M1,p(range(n))

def minimum_upset_random(M,N=10):
    """
    EXAMPLES:
        sage: M = matrix(QQ,[
        [0 , 0 , 1  , 0 , 0 , 0 ],
        [1,   0 ,  0,  1,  1,   0  ],
        [0 , 1 ,  0 ,  1 , 1  , 0  ],
        [1 , 0 , 0,  0 ,  0 , 0  ],
        [1 , 0 , 0 , 1 , 0 , 0  ],
        [1 ,  1  ,  1  , 1  , 1  , 0 ]
        ])
        sage: minimum_upset_random(M)
        (
        [0 0 1 1 0 1]                    
        [1 0 0 1 0 1]                    
        [0 1 0 0 0 0]                    
        [0 0 1 0 0 0]                    
        [1 1 1 1 0 1]                    
        [0 0 1 1 0 0], [1, 2, 0, 3, 5, 4]
        )

    """
    n = len(M.rows())
    Sn = SymmetricGroup(n)
    M1 = M
    wins = sum([sum([M1[j][i] for i in range(j,6)]) for j in range(6)])
    g0 = Sn(1)
    for k in range(N):
        g = Sn.random_element()
        P = g.matrix()
        M0 = P*M1*P^(-1)
        if sum([sum([M0[j][i] for i in range(j,6)]) for j in range(6)])>wins:
            M1 = M0
            g0 = g*g0
    return M1,g0(range(n))

def elo_rating(A):
    """
    A is a signed adjacency matrix for a directed graph.

    Returns elo ratings of the vertices of Gamma = Graph(A) 
        
    EXAMPLES:
        sage: A = matrix(QQ,[
        [0 , -1 , 1  , -1 , -1 , -1 ],
        [1,   0 ,  -1,  1,  1,   -1  ],
        [-1 , 1 ,  0 ,  1 , 1  , -1  ],
        [1 , -1 , -1,  0 ,  -1 , -1  ],
        [1 , - 1 , - 1 , 1 , 0 , - 1  ],
        [1 ,  1  ,  1  , 1  , 1  , 0 ]
        ])
        sage: elo_rating(A)
        (85.124, 104.79, 104.88, 85.032, 94.876, 124.53)

    """
    n = len(A.rows())
    RR = RealField(prec=20)
    V = RR^n
    K = 10
    r0 = 100 # initial rating
    r = n*[r0]
    for i in range(n):
        for j in range(n):
            if i<>j and A[i][j]==1:
                S = 1
            elif i<>j and A[i][j]==-1:
                S = 0
            else:
                S = 1/2
            mu = 1/(1+e^(-(r[i]-r[j])/400))
            r[i] = r[i]+K*(S-mu)
    return V(r)

def triangles(Gamma):
    """
    Returns a list of distinct (unoriented) triangles in Gamma,
    represented as 3-tuples of vertices.

    EXAMPLES:
        sage: A = matrix(QQ,[
        [0 , -1 , 1  , -1 , -1 , -1 ],
        [1,   0 ,  -1,  1,  1,   -1  ],
        [-1 , 1 ,  0 ,  1 , 1  , -1  ],
        [1 , -1 , -1,  0 ,  -1 , -1  ],
        [1 , - 1 , - 1 , 1 , 0 , - 1  ],
        [1 ,  1  ,  1  , 1  , 1  , 0 ]
        ])
        sage: Gamma = Graph(A,format="weighted_adjacency_matrix")
        sage: Gamma.triangles_count()                            
        20
        sage: T = triangles(Gamma)
        sage: len(T)
        20

    """
    E = Gamma.edges()
    V = Gamma.vertices()
    N = len(V)
    E0 = [(x[0],x[1]) for x in E]
    T = []
    for i1 in range(N):
        v1 = V[i1]
        for i2 in range(i1,N):
            v2 = V[i2]
            for i3 in range(i2,N):
                v3 = V[i3]
                #print v1,v2,v3
                if ((v1,v2) in E0) or ((v2,v1) in E0):
                    if ((v2,v3) in E0) or ((v3,v2) in E0):
                        if ((v1,v3) in E0) or ((v3,v1) in E0):
                            T.append([v1,v2,v3])
    return T

def grad(Gamma, f):
    """
    This returns the gradient of the function f : V -> ZZ, as a
    skew-symmetric matrix M with respect to the graph Gamma.

    EXAMPLE:
        sage: Gamma = graphs.DiamondGraph()
        sage: f = [1, 2, -3, -4]
        sage: grad(Gamma, f)
        [ 0 -1  4  5]
        [ 1  0  5  6]
        [-4 -5  0  1]
        [-5 -6 -1  0]

    """
    V = Gamma.vertices()
    n = len(V)
    M = matrix(ZZ, [[f[j]-f[i] for i in V] for j in V])
    return M


def curl(Gamma, M, i,j,k):
    """
    This returns the curl of the skew-symmetric matrix M with
    respect to the graph Gamma.

    EXAMPLE:
        sage: A = matrix(QQ,[
        [0 , -1 , 1  , -1 , -1 , -1 ],
        [1,   0 ,  -1,  1,  1,   -1  ],
        [-1 , 1 ,  0 ,  1 , 1  , -1  ],
        [1 , -1 , -1,  0 ,  -1 , -1  ],
        [1 , - 1 , - 1 , 1 , 0 , - 1  ],
        [1 ,  1  ,  1  , 1  , 1  , 0 ]
        ])
        sage: Gamma = Graph(A,format="weighted_adjacency_matrix")
        sage: curl(Gamma, A, 1, 4, 5)
        1
        sage: curl(Gamma, A, 1, 2, 4)
        -1
        sage: curl(Gamma, A, 0, 2, 3)
        3
        sage: Gamma = graphs.DiamondGraph()
        sage: f = [1, 2, -3, -4]
        sage: grad(Gamma, f)
        [ 0 -1  4  5]
        [ 1  0  5  6]
        [-4 -5  0  1]
        [-5 -6 -1  0]
        sage: curl(Gamma, grad(Gamma, f), 1,2,3)
        0

    """
    T = triangles(Gamma)
    if ([i,j,k] in T) or ([j,k,i] in T) or ([k,i,j] in T):
        return M[i][j]+M[j][k]+M[k][i]
    return 0

def div(Gamma, M, i):
    """
    This returns the divergence of the skew-symmetric matrix M with
    respect to the graph Gamma.

    EXAMPLE:
        sage: Gamma = graphs.DiamondGraph()
        sage: f = [1, 2, -3, -4]
        sage: M = grad(Gamma, f)
        sage: [div(Gamma, M, i) for i in range(4)]
        [8, 12, -8, -12]

    """
    n = len(M.rows())
    return sum([M[i][j] for j in range(n)])

def curl_matroid(M, A, i,j,k):
    """
    This returns the curl of the nxn skew-symmetric matrix A with
    respect to the matroid M, where n = size(M) (size of the groundset
    of M).

    EXAMPLE:
        sage: M0 = matrix(GF(2), [[1,0,0,0,1,1,1],[0,1,0,1,0,1,1],[0,0,1,1,1,0,1]])
        sage: M = sage.matroids.linear_matroid.LinearMatroid(M0)
        sage: CM = [x for x in M.circuits()]
        sage: TM = [x for x in CM if len(x)==3]
        sage: TM
        [frozenset({1, 2, 3}),
         frozenset({0, 2, 4}),
         frozenset({3, 4, 5}),
         frozenset({0, 1, 5}),
         frozenset({2, 5, 6}),
         frozenset({1, 4, 6}),
         frozenset({0, 3, 6})]
        sage: M.size()
        7
        sage: M77 = MatrixSpace(ZZ, 7, 7)
        sage: B = M77.random_element()
        sage: C = B - B.transpose()
        sage: C
        [  0   2  -1   0   0   1 -12]
        [ -2   0  -5   0   3  -1  10]
        [  1   5   0  -1   4   4  21]
        [  0   0   1   0   0  40  -1]
        [  0  -3  -4   0   0  -5  -3]
        [ -1   1  -4 -40   5   0  -1]
        [ 12 -10 -21   1   3   1   0]
        sage: curl_matroid(M, C, 1, 4, 6)
        -10
        sage: curl_matroid(M, C, 1, 2, 4)
        0
        sage: curl_matroid(M, C, 1, 2, 3)
        -6

    """
    CM = [x for x in M.circuits()]
    T = [x for x in CM if len(x)==3]
    if ({i,j,k} in T) or ({j,k,i} in T) or ({k,i,j} in T):
        return A[i][j]+A[j][k]+A[k][i]
    return 0

def linearly_independent_columns(A):
    """
    INPUT: A = mxn matrix.
    OUTPUT: The list of r-tuples of columns of independent
            columns, where r is the rank of A

    EXAMPLES:
        sage: A = matrix(GF(2), [[1,0,0,0,1,1,1],[0,1,0,1,0,1,1],[0,0,1,1,1,0,1]])
        sage: A
        [1 0 0 0 1 1 1]
        [0 1 0 1 0 1 1]
        [0 0 1 1 1 0 1]
        sage: linearly_independent_columns(A)
        [[0, 1, 2],
         [0, 1, 3],
         [0, 1, 4],
         [0, 1, 6],
         [0, 2, 3],
         [0, 2, 5],
         [0, 2, 6],
         [0, 3, 4],
         [0, 3, 5],
         [0, 4, 5],
         [0, 4, 6],
         [0, 5, 6],
         [1, 2, 4],
         [1, 2, 5],
         [1, 2, 6],
         [1, 3, 4],
         [1, 3, 5],
         [1, 3, 6],
         [1, 4, 5],
         [1, 5, 6],
         [2, 3, 4],
         [2, 3, 5],
         [2, 3, 6],
         [2, 4, 5],
         [2, 4, 6],
         [3, 4, 6],
         [3, 5, 6],
         [4, 5, 6]]


    """
    L = []
    colsA = A.columns()
    r = A.rank()
    n = len(colsA)
    Cnr = Combinations(n,r)
    for c in Cnr:
        B = matrix([colsA[i] for i in c])
	if r==B.rank():
	    L=L+[c]
    return L

def codim1_sublists(L):
    """
    L is a list of tuples, all of the same size. 
    This function returns the sublists of one smaller size.

    EXAMPLES:
        sage: L = [[0, 1, 2], [0, 1, 3]]
        sage: codim1_sublists(L)
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
    """
    n = len(L[0])
    L0 = []
    for x in L:
        L0=L0+Combinations(x,n-1).list()
    L1 = []
    #print L0,uniq(L0)
    for x in L0:
        L1 = L1+[tuple(x)]
    return uniq(L1)

def subfaces_graph(L):
    """
    L is a list of facets.
    This fcn constructs the bipartite graph whose
    vertices are L and its ridges R. Vertices v,w are
    connected by an edge iff v in L and w in R and
    w is incident to v (ie, contained in v as a set).

    EXAMPLES:
        sage: L = [[0, 1, 2], [0, 1, 3]]
        sage: codim1_sublists(L)
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]
        sage: subfaces_graph(L)
        Graph on 7 vertices
        sage: L = [[0, 1, 2], [0, 1, 3], [0, 1, 4], [0, 1, 6], [0, 2, 3], [0, 2, 5 ], [0, 2, 6], [0, 3, 4], [0, 3, 5 ], [0, 4, 5 ], [0, 4, 6], [0, 5, 6], [1, 2, 4], [1, 2, 5 ], [1, 2, 6], [1, 3, 4], [1, 3, 5 ], [1, 3, 6], [1, 4, 5 ], [1, 5, 6], [2, 3, 4], [2, 3, 5 ], [2, 3, 6], [2, 4, 5 ], [2, 4, 6], [3, 4, 6], [3, 5, 6], [4, 5, 6]] # max ind sets for Fano matroid
        sage: Gamma = subfaces_graph(L)
        sage: V = Gamma.vertices()
        sage: pos_dict = {}
        sage: for i in range(16):
            x = float(10*cos(pi/2 + ((2*pi)/16)*i))
            y = float(10*sin(pi/2 + ((2*pi)/16)*i))
            pos_dict[V[i]] = [x,y]
        ....:     
        sage: Gamma.show(pos=pos_dict)
        sage: pos_dict = {}
        sage: for i in range(21):
            x = float(1000*cos(pi/2 + ((2*pi)/21)*i))
            y = float(1000*sin(pi/2 + ((2*pi)/21)*i))
            pos_dict[Vsort[i]] = [x,y]
        ....:     
        sage: for i in range(21,35):
            x = float(500*cos(pi/2 + ((2*pi)/14)*i))
            y = float(500*sin(pi/2 + ((2*pi)/14)*i))
            pos_dict[Vsort[i]] = [x,y]
        ....:     
        sage: for i in range(35,49):
            x = float(1500*cos(pi/2 + ((2*pi)/14)*i))
            y = float(1500*sin(pi/2 + ((2*pi)/14)*i))
            pos_dict[Vsort[i]] = [x,y]
        ....:     
        Gamma.show(pos=pos_dict,vertex_size=500, dpi=400)
    """
    L0 = [tuple(x) for x in L]
    L1 = codim1_sublists(L0)
    V = L0+L1
    E = []
    for v in L0:
        for w in L1:
	    if Set(w).issubset(Set(v)):
	        E = E+[(v,w)]
    return Graph([V,E])
    

def simplex_laplacian_combinatorial(X,i):
    """
    Returns the i-th combinatorial Laplacian of the simplicial complex X.

    EXAMPLES:
        sage: X = simplicial_complexes.Sphere(3)
        sage: simplex_laplacian_combinatorial(X,1)
        [5 0 0 0 0 0 0 0 0 0]
        [0 5 0 0 0 0 0 0 0 0]
        [0 0 5 0 0 0 0 0 0 0]
        [0 0 0 5 0 0 0 0 0 0]
        [0 0 0 0 5 0 0 0 0 0]
        [0 0 0 0 0 5 0 0 0 0]
        [0 0 0 0 0 0 5 0 0 0]
        [0 0 0 0 0 0 0 5 0 0]
        [0 0 0 0 0 0 0 0 5 0]
        [0 0 0 0 0 0 0 0 0 5]

    """
    C = X.chain_complex(dimensions=[i-1,i,i+1])    
    return C.differential(i+1)*C.differential(i+1).transpose()+C.differential(i).transpose()*C.differential(i)

def simplex_laplacian_combinatorial_up(X,i):
    """
    Returns the i-th combinatorial Laplacian of the simplicial complex X.

    EXAMPLES:
        sage: X = simplicial_complexes.Sphere(3)
        sage: simplex_laplacian_combinatorial_up(X,1)
        [ 3  1 -1 -1 -1 -1  1  0  0  0]
        [ 1  3  0 -1  0  1 -1 -1  0  1]
        [-1  0  3 -1 -1  0 -1  0 -1 -1]
        [-1 -1 -1  3 -1  0  0 -1  0  1]
        [-1  0 -1 -1  3  1  0  1  1  0]
        [-1  1  0  0  1  3  1 -1 -1  0]
        [ 1 -1 -1  0  0  1  3  0 -1 -1]
        [ 0 -1  0 -1  1 -1  0  3 -1  1]
        [ 0  0 -1  0  1 -1 -1 -1  3 -1]
        [ 0  1 -1  1  0  0 -1  1 -1  3]
        sage: S = SimplicialComplex(maximal_faces=[(1,2,3), (1,2,4),\
         (1,2,5), (1,3,4),(1,3,5),(2,3,4),(2,3,5)])
        sage: simplex_laplacian_combinatorial_up(S,1)
        [ 3 -1 -1 -1  1  1  1  0  0]
        [-1  3 -1 -1 -1  0  0  1  1]
        [-1 -1  2  0  0 -1  0 -1  0]
        [-1 -1  0  2  0  0 -1  0 -1]
        [ 1 -1  0  0  3 -1 -1  1  1]
        [ 1  0 -1  0 -1  2  0 -1  0]
        [ 1  0  0 -1 -1  0  2  0 -1]
        [ 0  1 -1  0  1 -1  0  2  0]
        [ 0  1  0 -1  1  0 -1  0  2]

    """
    C = X.chain_complex(dimensions=[i,i+1])
    Q = C.differential(i+1)*C.differential(i+1).transpose()
    return Q

def simplex_degree_configuration(X, i):
    """
    Returns the i-th degree configuration of the simplicial complex X.

    EXAMPLES:
        sage: X = simplicial_complexes.Sphere(3)
        sage: [2, 2, 2, 2, 2, 2]

    """
    Q = simplex_laplacian_combinatorial_up(X,i)
    d = Q.diagonal()
    return d
    
def simplex_laplacian_combinatorial_down(X,i):
    """
    Returns the i-th combinatorial Laplacian of the simplicial complex X.

    EXAMPLES:
        sage: X = simplicial_complexes.Sphere(3)
        sage: simplex_laplacian_combinatorial_down(X,1)
        [ 2 -1  1  1  1  1 -1  0  0  0]
        [-1  2  0  1  0 -1  1  1  0 -1]
        [ 1  0  2  1  1  0  1  0  1  1]
        [ 1  1  1  2  1  0  0  1  0 -1]
        [ 1  0  1  1  2 -1  0 -1 -1  0]
        [ 1 -1  0  0 -1  2 -1  1  1  0]
        [-1  1  1  0  0 -1  2  0  1  1]
        [ 0  1  0  1 -1  1  0  2  1 -1]
        [ 0  0  1  0 -1  1  1  1  2  1]
        [ 0 -1  1 -1  0  0  1 -1  1  2]

    """
    C = X.chain_complex(dimensions=[i-1,i,i+1])
    Q = C.differential(i).transpose()*C.differential(i)
    return Q

