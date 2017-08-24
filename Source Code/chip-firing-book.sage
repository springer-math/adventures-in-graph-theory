"""
Implements some aspects of chip-firing game functions.

Thanks to David Perkinson, at Reed College, Sage now has excellent 
chip-firing functions, see
http://doc.sagemath.org/html/en/thematic_tutorials/sandpile.html

The functions in this module were written independently.
Used in [1].

REFERENCES:
    [0] N. Biggs, "Chip firing and the critical group of a 
    graph," J. Algeb. Comb. 9 (1999)25-45.
    [1] D. Joyner and C. Melles, "Adventures in graph theory,"
    Birkhauser, to appear.

Created, 2010-04. Copyright David Joyner, 2017-06-08
license: modified BSD
"""

def active_vertices(G, s):
    """
    Returns the list of active vertices.

    INPUT:
     G - a graph
     s - a configuration (implemented as a list 
                          or a dictionary keyed on
                          the vertices of the graph)

    EXAMPLES:
        sage: A = matrix([[0,1,1,0,0],[1,0,1,0,0],[1,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]])
        sage: G = Graph(A, format = "adjacency_matrix", weighted = True)
        sage: s = {0: 3, 1: 1, 2: 0, 3: 1, 4: 1}
        sage: active_vertices(G, s)
        [0, 4]

    """
    V = G.vertices()
    degs = [G.degree(v) for v in V]
    active = [v for v in V if degs[V.index(v)]<=s[v]]
    return active

def active_sets(G, s, source = 0):
    """
    Returns the list of active vertices.

    INPUT:
     G - a graph
     s - a configuration (implemented as a list 
                          or a dictionary keyed on
                          the vertices of the graph)

    EXAMPLES:
        sage: A = matrix([[0,1,1,0,0],[1,0,1,0,0],[1,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]])
        sage: G = Graph(A, format = "adjacency_matrix", weighted = True)
        sage: s = {0: 3, 1: 1, 2: 0, 3: 1, 4: 1}
        sage: active_vertices(G, s)
        [0, 4]

    """
    V = G.vertices()
    n = len(V)
    IS = IntegerModRing(2)^(n-1)
    ac = []
    vs = vector(s.values())
    Q = G.laplacian_matrix()
    for x in IS:
        vx = vector([0]+[int(z) for z in x]) # assumes source = 0
        y = vs - Q*vx
        add_ac = True
        #print y
        for i in range(len(y.list())):
            if i<>source and y[i]<0:
                add_ac = False
        if add_ac == True:
            ac.append(vx.list())
    return ac

def stable_vertices(G, s, source = None):
    """
    Returns the list of stable vertices.

    INPUT:
     G - a graph
     s - a configuration (implemented as a list 
                          or a dictionary keyed on
                          the vertices of the graph)

    EXAMPLES:
        sage: A = matrix([[0,1,1,0,0],[1,0,1,0,0],[1,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]])
        sage: G = Graph(A, format = "adjacency_matrix", weighted = True)
        sage: s = {0: 3, 1: 1, 2: 0, 3: 1, 4: 1}
        sage: stable_vertices(G, s)
        [1, 2, 3]

    """
    V = G.vertices()
    degs = [G.degree(v) for v in V]
    if source==None:
        stable = [v for v in V if degs[V.index(v)]>s[v]]
    else:
        stable = [v for v in V if degs[V.index(v)]>s[v] and v!=source]
    return stable

def fire(G, s, v0):
    """
    Returns the configuration after firing the active vertex v.

    INPUT:
     G - a graph
     s - a configuration (implemented as a list 
                          or a dictionary keyed on
                          the vertices of the graph)
     v - a vertex of the graph

    EXAMPLES:
        sage: A = matrix([[0,1,1,0,0],[1,0,1,0,0],[1,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]])
        sage: G = Graph(A, format = "adjacency_matrix", weighted = True)
        sage: s1 = {0: 3, 1: 1, 2: 0, 3: 1, 4: 1}
        sage: s2 = fire(G, s1, 0); s2
        {0: 1, 1: 2, 2: 1, 3: 1, 4: 1}
        sage: s3 = fire(G, s2, 1); s3
        {0: 2, 1: 0, 2: 2, 3: 1, 4: 1}
        sage: s4 = fire(G, s3, 4); s4
        {0: 2, 1: 0, 2: 2, 3: 2, 4: 0}
        sage: s5 = fire(G, s4, 3); s5
        {0: 2, 1: 0, 2: 3, 3: 0, 4: 1}
        sage: s6 = fire(G, s5, 2); s6
        {0: 3, 1: 1, 2: 0, 3: 1, 4: 1}

    """
    V = G.vertices()
    j = V.index(v0)
    s1 = copy(s)
    if not(v0 in V):
        raise ValueError, "the last argument must be a vertex of the graph."
    if not(v0 in active_vertices(G, s)):
        raise ValueError, "the last argument must be an active vertex of the graph."
    degs = [G.degree(w) for w in V]
    for w in V:
         if w == v0:
             s1[v0] = s[v0] - degs[j]
         if w in G.neighbors(v0):
             s1[w] = s[w]+1
    return s1

def fire_set(G, s, S0):
    """
    Returns the configuration after firing the active vertex v.

    INPUT:
     G - a graph
     s - a configuration (implemented as a list 
                          or a dictionary keyed on
                          the vertices of the graph)
     S - a set of vertices of the graph

    EXAMPLES:


    """
    V = G.vertices()
    Q = G.laplacian_matrix()
    n = len(V)
    s1 = copy(s)
    for v0 in S0:
        if not(v0 in V):
            raise ValueError, "the last argument must be a vertex of the graph."
    vs = vector(s.values())
    ss = vs - Q*vector(S0)
    #print ss
    Ls = ss.list()
    s1 = dict([(i,Ls[i]) for i in range(n)])
    return s1


def stabilize(G, s, source, legal_sequence = False):
    """
    Returns the stable configuration of the graph originating from
    the given legal configuration. If legal_sequence = True then the
    sequence of firings is also returned. By van den Heuvel [1],
    the number of firings needed to compute a critical configuration 
    is < 3(S+2|E|)|V|^2, where S is the sum of the positive 
    weights in the configuration.

    EXAMPLES:
        sage: A = matrix([[0,1,1,0,0],[1,0,1,0,0],[1,1,0,1,0],[0,0,1,0,1],[0,0,0,1,0]])
        sage: G = Graph(A, format="weighted_adjacency_matrix")
        sage: s = {0: 3, 1: 1, 2: 0, 3: 1, 4: -5}
        sage: stabilize(G, s, 4)
        {0: 0, 1: 1, 2: 2, 3: 1, 4: -4}
        sage: stabilize(G, s, 4,legal_sequence = True)
        ({0: 0, 1: 1, 2: 2, 3: 1, 4: -4}, [0, 1, 0, 2, 1, 0, 3, 2, 1, 0])
        sage: Gamma = graphs.CycleGraph(4)
        sage: Gamma.add_edge(1,3)
        sage: critical_group(Gamma)
        Finitely generated module V/W over Integer Ring with invariants (8)
        sage: s = {0: -2, 1: 1, 2: 1, 3: 0}
        sage: stabilize(Gamma, s, 0) == s
        True
        sage: s = {0: -2, 1: 0, 2: 1, 3: 1}
        sage: stabilize(Gamma, s, 0) == s
        True
        sage: s = {0: -3, 1: 1, 2: 1, 3: 1}
        sage: stabilize(Gamma, s, 0) == s
        True
        sage: s = {0: -3, 1: 1, 2: 0, 3: 2}
        sage: stabilize(Gamma, s, 0) == s
        True
        sage: s = {0: -3, 1: 2, 2: 0, 3: 1}
        sage: stabilize(Gamma, s, 0) == s
        True
        sage: s = {0: -4, 1: 2, 2: 0, 3: 2}
        sage: stabilize(Gamma, s, 0) == s
        True
        sage: s = {0: -4, 1: 1, 2: 1, 3: 2}
        sage: stabilize(Gamma, s, 0) == s
        True
        sage: s = {0: -4, 1: 2, 2: 1, 3: 1}
        sage: stabilize(Gamma, s, 0) == s
        True
        sage: s = {0: -5, 1: 2, 2: 1, 3: 2}
        sage: stabilize(Gamma, s, 0) == s
        True

    REFERENCES:
        [1] J. van den Heuvel, "Algorithmic aspects of a chip-firing
            game," preprint.
    """
    V = G.vertices()
    E = G.edges()
    fire_number = 3*len(V)^2*(sum([s[v] for v in V if s[v]>0])+2*len(E))+len(V)
    if legal_sequence:
        seq = []
    stab = []
    ac = active_vertices(G,s)
    for i in range(fire_number):
        if len(ac)>0:
            s = fire_set(G,s,ac)
            if legal_sequence:
                seq.append(ac)
        else:
            stab.append(s)
            break
        ac = active_vertices(G,s)
    if len(stab)==0:
        raise ValueError, "No stable configuration found."
    if legal_sequence:
        return stab[0], seq
    else:
        return stab[0]

def stabilize2(G, s, source, legal_sequence = False):
    """
    Returns the stable configuration of the graph originating from
    the given legal configuration. If legal_sequence = True then the
    sequence of firings is also returned. By van den Heuvel [1],
    the number of firings needed to compute a critical configuration 
    is < 3(S+2|E|)|V|^2, where S is the sum of the positive 
    weights in the configuration.

    EXAMPLES:
        sage: Gamma = graphs.WheelGraph(4)
        sage: s = {0: -5, 1: 2, 2: 1, 3: 2}
        sage: stabilize2(Gamma, s, 0)
        {0: -2, 1: 1, 2: 0, 3: 1}


    REFERENCES:
        [1] J. van den Heuvel, "Algorithmic aspects of a chip-firing
            game," preprint.
    """
    V = G.vertices()
    E = G.edges()
    fire_number = 3*len(V)^2*(sum([s[v] for v in V if s[v]>0])+2*len(E))+len(V)
    if legal_sequence:
        seq = []
    stab = []
    for j in range(fire_number):
        ac = active_sets(G,s,source=source)
        ac.sort()
        ac.reverse()
        #print ac
        if len(ac)>0:
            S = ac[0]
            if vector(S)<>vector(S)*0:
                #print i, s
                if s <> fire_set(G,s,S):
                    s = fire_set(G,s,S)
                else:
                    break
                ac = active_sets(G,s,source=source)
                #print i, s
                if legal_sequence:
                    seq.append(S)
    stab.append(s)
    if len(stab)==0:
        raise ValueError, "No stable configuration found."
    if legal_sequence:
        return stab[0], seq
    else:
        return stab[0]

def add_configurations(G, s1, s2):
    """
    Returns the critical configuration associated to the sum of s1, s2.

    INPUT:
     G - a graph
     s1, s2 - two critcal configuration (implemented as a list 
                          or a dictionary keyed on
                          the vertices of the graph)

    EXAMPLES:
       sage: Gamma = graphs.WheelGraph(6)
       sage: s = {0: -4, 1: 0, 2: 2, 3: 0, 4: 1, 5: 1}
       sage: add_configurations(Gamma, s, s)
       {0: -7, 1: 1, 2: 1, 3: 1, 4: 2, 5: 2}
       sage: add_configurations(Gamma, s, s).values()
       [-7, 1, 1, 1, 2, 2]
       sage: s2 = add_configurations(Gamma, s, s); s2.values()
       [-7, 1, 1, 1, 2, 2]
       sage: s3 = add_configurations(Gamma, s2, s); s3.values()
       [-6, 0, 2, 0, 2, 2]
       sage: s4 = add_configurations(Gamma, s3, s); s4.values()
       [-7, 2, 1, 2, 1, 1]
       sage: s5 = add_configurations(Gamma, s4, s); s5.values()
       [-6, 1, 2, 1, 1, 1]
       sage: s6 = add_configurations(Gamma, s5, s); s6.values()
       [-9, 2, 1, 2, 2, 2]
       sage: s7 = add_configurations(Gamma, s6, s); s7.values()
       [-8, 1, 2, 1, 2, 2]
       sage: s8 = add_configurations(Gamma, s7, s); s8.values()
       [-6, 1, 0, 1, 2, 2]
       sage: s9 = add_configurations(Gamma, s8, s); s9.values()
       [-8, 2, 2, 2, 1, 1]
       sage: s10 = add_configurations(Gamma, s9, s); s10.values()
       [-6, 2, 0, 2, 1, 1]
       sage: s11 = add_configurations(Gamma, s10, s); s11.values()
       [-10, 2, 2, 2, 2, 2]

    """
    v1 = vector(ZZ, s1.values())
    v2 = vector(ZZ, s2.values())
    v = v1+v2
    n = len(s2.values())
    s = {x: list(v)[x] for x in range(n)}
    return stabilize(G, s, 0)

def laplacian_moore_penrose_inverse(Gamma):
    """
    Returns the Moore-Penrose inverse of the Laplacian,
    using Shokreih's formula.

    INPUT: Graph, Gamma
    OUTPUT: Q^+

    EXAMPLES:
        sage: Gamma = graphs.CycleGraph(4)
        sage: Gamma.add_edge(1,3)
        sage: laplacian_moore_penrose_inverse(Gamma)
        [ 5/16 -1/16 -3/16 -1/16]
        [-1/16  3/16 -1/16 -1/16]
        [-3/16 -1/16  5/16 -1/16]
        [-1/16 -1/16 -1/16  3/16]

    """
    Q = Gamma.laplacian_matrix()
    V = Gamma.vertices()
    n = len(V)
    MS = MatrixSpace(ZZ, n, n)
    J = MS(n*n*[1])
    Qplus = (Q+(1/n)*J)^(-1)-(1/n)*J
    return Qplus

def shokreih_pairing(Gamma, s1, s2):
    """
    Returns the monodromy pairing of two critical configurations,
    using Shokreih's formula.

    INPUT: critical configurations, s1, s2
    OUTPUT: [s1, s2] (in QQ).

    EXAMPLES:
        sage: Gamma = graphs.CycleGraph(4)
        sage: Gamma.add_edge(1,3)
        sage: Qplus = laplacian_moore_penrose_inverse(Gamma)
        sage: c0 = vector([-5,2,1,2])
        sage: c1 = vector([-4,2,1,1])
        sage: c2 = vector([-4,1,1,2])
        sage: c3 = vector([-4,2,0,2])
        sage: c4 = vector([-3,2,0,1])
        sage: c5 = vector(ZZ, [-3, 1, 0, 2])
        sage: c6 = vector([-3,0,1,2])
        sage: c7 = vector([-3,2,1,0])
        sage: c2*Qplus*c5
        49/8
        sage: SP = [[shokreih_pairing(Gamma, cc[i], cc[j]) for i in range(8)] for j in range(8)]
        sage: matrix(SP)
        [  0 1/2 1/2   0 1/2 1/2   0   0]
        [1/2 5/8 3/8   0 1/8 7/8 1/4 3/4]
        [1/2 3/8 5/8   0 7/8 1/8 3/4 1/4]
        [  0   0   0   0   0   0   0   0]
        [1/2 1/8 7/8   0 5/8 3/8 1/4 3/4]
        [1/2 7/8 1/8   0 3/8 5/8 3/4 1/4]
        [  0 1/4 3/4   0 1/4 3/4 1/2 1/2]
        [  0 3/4 1/4   0 3/4 1/4 1/2 1/2]

    """
    Q = Gamma.laplacian_matrix()
    V = Gamma.vertices()
    n = len(V)
    MS = MatrixSpace(ZZ, n, n)
    J = MS(n*n*[1])
    Qplus = laplacian_moore_penrose_inverse(Gamma)
    a = s1*Qplus*s2
    return a-floor(a)

def reduced_configurations0(Gamma, weight = 0):
    """
    Returns the list of all reduced configurations of Gamma.


    EXAMPLES:
        sage: Gamma = graphs.WheelGraph(4)
        sage: L = reduced_configurations(Gamma)
        sage: len(L)
        16
        sage: critical_group(Gamma)
        Finitely generated module V/W over Integer Ring with invariants (4, 4)
        sage: Gamma = graphs.WheelGraph(5)
        sage: len(reduced_configurations(Gamma))
        45
        sage: critical_group(Gamma)
        Finitely generated module V/W over Integer Ring with invariants (3, 15)

    """
    degs = Gamma.degree_sequence()
    d = max(degs)
    n = len(degs)
    C = IntegerModRing(int(d))^(n-1)
    L = []
    for c in C:
        if sum(c)==weight:
            c = [int(x) for x in c]
            s = [-sum(c)]+c
            sd = [(i,s[i]) for i in range(n)]
            sr = stabilize2(Gamma, dict(sd), source=0, legal_sequence = False)
            #print c, sd, sr
            L.append(sr.values())
    Is = [L.index(x) for x in L]
    Is = Set(Is)
    LL = [L[i] for i in Is]
    LL.sort()
    return LL


def reduced_configurations(Gamma):
    """
    Returns the list of all reduced configurations of Gamma.


    EXAMPLES:
        sage: Gamma = graphs.WheelGraph(4)
        sage: L = reduced_configurations(Gamma)
        sage: len(L)
        16
        sage: critical_group(Gamma)
        Finitely generated module V/W over Integer Ring with invariants (4, 4)
        sage: Gamma = graphs.WheelGraph(5)
        sage: len(reduced_configurations(Gamma))
        45
        sage: critical_group(Gamma)
        Finitely generated module V/W over Integer Ring with invariants (3, 15)

    """
    degs = Gamma.degree_sequence()
    d = max(degs)
    n = len(degs)
    C = IntegerModRing(int(d))^(n-1)
    L = []
    for c in C:
        c = [int(x) for x in c]
        s = [-sum(c)]+c
        sd = [(i,s[i]) for i in range(n)]
        sr = stabilize2(Gamma, dict(sd), source=0, legal_sequence = False)
        #print c, sd, sr
        L.append(sr.values())
    Is = [L.index(x) for x in L]
    Is = Set(Is)
    LL = [L[i] for i in Is]
    LL.sort()
    return LL

def linear_system(D, Gamma):
    """
    Returns linear system attached to the divisor D.

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: is_harmonic_graph_morphism(Gamma1, Gamma2, f)
        True
        sage: PhiV = matrix_of_graph_morphism_vertices(Gamma1, Gamma2, f); PhiV
        [1 0 1 0]
        [0 1 0 1]
        sage: D = vector([1,0,0,1])
        sage: PhiV*D
        (1, 1)
        sage: linear_system(PhiV*D, Gamma1)
        [(2, 0), (1, 1), (0, 2)]
        sage: linear_system(D, Gamma2)
        [(0, 2, 0, 0), (0, 0, 2, 0), (1, 0, 0, 1)]
        sage: [PhiV*x for x in linear_system(D, Gamma2)]
        [(0, 2), (2, 0), (1, 1)]

    """
    Q = Gamma.laplacian_matrix()
    CS = Q.column_space()
    N = len(D.list())
    d = sum(D.list())
    #print d
    lin_sys = []
    if d < 0:
        return lin_sys
    if (d == 0) and (D in CS):
        lin_sys = [CS(0)]
        return lin_sys
    elif (d == 0):
        return lin_sys
    S = IntegerModRing(d+1)^N
    V = QQ^N
    for v in S:
        v = V(v)
        #print D-v,v,D
        if D-v in CS:
            lin_sys.append(v)
    return lin_sys

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

def effective_nonspecial_divisors(Gamma):
    """
    A divisor D on a graph G is called non-special if 
        1. deg(D)=g-1,where g is the genus of G, and
        2. linear_system(D) = [], i.e., the linear system associated to D is empty.
    ### This is a silly function ... ###

    EXAMPLES:
        sage: Gamma2 = graphs.WheelGraph(4)
        sage: genus(Gamma2)
        3
        sage: effective_nonspecial_divisors(Gamma2)
        []
        sage: Gamma2 = graphs.WheelGraph(5)
        sage: genus(Gamma2)
        4
        sage: len(effective_nonspecial_divisors(Gamma2))
        0
        sage: effective_nonspecial_divisors(Gamma2)
        []
    """
    ns_divs = []
    d = genus(Gamma)-1
    N = len(Gamma.vertices())
    S = IntegerModRing(d+1)^N
    V = QQ^N
    for v in S:
        v = V(v)
        if sum(v.list())==d and linear_system(v, Gamma)==[]:
            ns_divs.append(v)
    return ns_divs

def canonical_divisor(Gamma):
    """
    EXAMPLES:
        sage: Gamma2 = graphs.WheelGraph(5)
        sage: canonical_divisor(Gamma2)
        (2, 1, 1, 1, 1)

    """
    V = Gamma.vertices()
    return vector([Gamma.degree(v)-2 for v in V])

def riemann_roch_dimension(D, Gamma):
    """
    The dimension r(D) of a linear system |D| is defined 
    as follows. 
    If |D| is empty, r(D) = -1. 
    Otherwise, if |D| is non-empty, r(D) is the largest 
    nonnegative integer d such that |D-D'| is non-empty, 
    for all effective divisors D' of degree d. 
    The below assumes deg(D)>0 and r(D)<=deg(D)

    EXAMPLES:
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: D = vector([1,0,0,1])
        sage: riemann_roch_dimension(D, Gamma2)
        1
        sage: riemann_roch_dimension(canonical_divisor(Gamma2), Gamma2)
        0
        sage: Gamma1 = graphs.TetrahedralGraph()
        sage: D = vector([1,0,0,1])
        sage: riemann_roch_dimension(D, Gamma1)
        0
        sage: D = vector([2,0,0,1])
        sage: riemann_roch_dimension(D, Gamma1)
        0
        sage: D = vector([3,0,0,1])
        sage: riemann_roch_dimension(D, Gamma1)
        1
        sage: Gamma2 = graphs.CubeGraph(2)
        sage: Gamma1 = Gamma2.subgraph(vertices = ['00', '01'], edges = [('00', '01')])
        sage: f = [['00', '01', '10', '11'], ['00', '01', '00', '01']]
        sage: D = vector([1,0,0,1])
        sage: PhiV = matrix_of_graph_morphism_vertices(Gamma1, Gamma2, f); PhiV
        [1 0 1 0]
        [0 1 0 1]
        sage: PhiV*D
        (1, 1)
        sage: linear_system(PhiV*D, Gamma1)
        [(2, 0), (1, 1), (0, 2)]
        sage: linear_system(D, Gamma2)
        [(0, 2, 0, 0), (0, 0, 2, 0), (1, 0, 0, 1)]
        sage: riemann_roch_dimension(D, Gamma2)
        1
        sage: riemann_roch_dimension(PhiV*D, Gamma1)
        2
        sage: genus(Gamma2)
        1
        sage: genus(Gamma1)
        0
        sage: D = vector([2,0,0,1])
        sage: riemann_roch_dimension(D, Gamma2)
        2
        sage: riemann_roch_dimension(PhiV*D, Gamma1)
        3
        sage: D = vector([2,1,0,1])
        sage: riemann_roch_dimension(D, Gamma2)
        3
        sage: riemann_roch_dimension(PhiV*D, Gamma1)
        4

        sage: Gamma2 = graphs.WheelGraph(4)
        sage: Gamma1 = graphs.WheelGraph(2)
        sage: genus(Gamma2)
        3
        sage: genus(Gamma1)
        0
        sage: f = [[0,1,2,3],[0,1,1,1]]
        sage: PhiV = matrix_of_graph_morphism_vertices(Gamma1, Gamma2, f); PhiV
        [1 0 0 0]
        [0 1 1 1]
        sage: D = vector([1,0,0,1])
        sage: riemann_roch_dimension(D, Gamma2)
        0
        sage: riemann_roch_dimension(PhiV*D, Gamma1)
        2
        sage: D = vector([2,0,0,1])
        sage: riemann_roch_dimension(D, Gamma2)
        0
        sage: riemann_roch_dimension(PhiV*D, Gamma1)
        3
        sage: D = vector([3,0,0,1])
        sage: riemann_roch_dimension(D, Gamma2)
        1
        sage: riemann_roch_dimension(PhiV*D, Gamma1)
        4
        sage: D = vector([2,1,1,-1])
        sage: riemann_roch_dimension(D, Gamma2)
        0
        sage: riemann_roch_dimension(PhiV*D, Gamma1)
        3
        sage: D = vector([3,1,1,-1])
        sage: riemann_roch_dimension(D, Gamma2)
        1
        sage: riemann_roch_dimension(PhiV*D, Gamma1)
        4
    """
    if linear_system(D, Gamma)==[]:
        return -1
    d = sum(D.list())
    N = len(Gamma.vertices())
    V = QQ^N
    for rD in range(d+1):
        S = IntegerVectors(rD+1,N)
        #S = IntegerModRing(rD+2)^N
        for v in S:
            v = V(v)
            if linear_system(D-v, Gamma)==[]:
            #if sum(v.list())==rD+1 and linear_system(D-v, Gamma)==[]:
                return rD
    return rD,"finished"

"""
Find small vector repns 

R1 = GF(5)^5
R2 = ZZ^5
for c1 in R1:
    c = R2(c1)
    v = c[0]*b0+c[1]*b1+c[2]*b2+c[3]*b3+c[4]*b4
    if is_stable(S, i, v, [2,3,6,7]):
        print c,v
#(2, 3, 0, 2, 2) (2, 3, 0, -5, 2, 2, -2, -2, 7)
#(3, 3, 0, 2, 2) (3, 3, 0, -6, 2, 2, -1, -2, 7)
#(4, 3, 0, 2, 2) (4, 3, 0, -7, 2, 2, 0, -2, 7)
#(2, 4, 0, 2, 2) (2, 4, 0, -6, 2, 2, -2, -2, 8)
#(3, 4, 0, 2, 2) (3, 4, 0, -7, 2, 2, -1, -2, 8)
#(4, 4, 0, 2, 2) (4, 4, 0, -8, 2, 2, 0, -2, 8)
#(2, 3, 1, 2, 2) (2, 3, 1, -6, 2, 2, -2, -3, 8)

for i1 in range(2):
  for i2 in range(2):
    for i3 in range(2):
      for i4 in range(2):
        for i5 in range(2):
          v = i1*BZ1[0]+i2*BZ1[1]+i3*BZ1[2]+i4*BZ1[3]+i5*BZ1[4]
          if not(v in B1):
                print i1,i2,i3,i4,i5,v

0 0 0 0 1 (0, 0, 0,  0, 0, 1, -1, -1, 1)
0 0 0 1 0 (0, 0, 0,  0, 1, 0, -1,  0, 1)
0 0 1 0 0 (0, 0, 1, -1, 0, 0,  0, -1, 1)
0 0 1 1 0 (0, 0, 1, -1, 1, 0, -1, -1, 2)
0 0 1 1 1 (0, 0, 1, -1, 1, 1, -2, -2, 3)
0 1 0 0 0 (0, 1, 0, -1, 0, 0,  0,  0, 1)
0 1 0 0 1 (0, 1, 0, -1, 0, 1, -1, -1, 2)
0 1 0 1 1 (0, 1, 0, -1, 1, 1, -2, -1, 3)
0 1 1 0 0 (0, 1, 1, -2, 0, 0,  0, -1, 2)
0 1 1 0 1 (0, 1, 1, -2, 0, 1, -1, -2, 3)
0 1 1 1 0 (0, 1, 1, -2, 1, 0, -1, -1, 3)
1 0 0 0 0 (1, 0, 0, -1, 0, 0,  1,  0, 0)
1 0 0 0 1 (1, 0, 0, -1, 0, 1,  0, -1, 1)
1 0 0 1 0 (1, 0, 0, -1, 1, 0,  0,  0, 1)
1 0 0 1 1 (1, 0, 0, -1, 1, 1, -1, -1, 2)
1 0 1 0 0 (1, 0, 1, -2, 0, 0,  1, -1, 1)
1 0 1 0 1 (1, 0, 1, -2, 0, 1,  0, -2, 2)
1 0 1 1 0 (1, 0, 1, -2, 1, 0,  0, -1, 2)
1 0 1 1 1 (1, 0, 1, -2, 1, 1, -1, -2, 3)
1 1 0 0 1 (1, 1, 0, -2, 0, 1,  0, -1, 2)
1 1 0 1 0 (1, 1, 0, -2, 1, 0,  0,  0, 2)
1 1 0 1 1 (1, 1, 0, -2, 1, 1, -1, -1, 3)
1 1 1 0 0 (1, 1, 1, -3, 0, 0,  1, -1, 2)
1 1 1 1 0 (1, 1, 1, -3, 1, 0,  0, -1, 3)
1 1 1 1 1 (1, 1, 1, -3, 1, 1, -1, -2, 4)

# via LLL: [ 0  0 -1  0  2  0  1 -1  1]

"""

def is_active(X, i, c, S0):
    """
    Returns True if the configuration c is stable,
    False otherwise.

     INPUT:
     X - a simplicial complex
     c - a configuration (implemented as a list having the same length 
          as the number of cols of the i-th boundary matrix of X)
     S0 - a set of source $i$-faces of X (ie, a spanning tree)


    EXAMPLES:
        sage: S = SimplicialComplex(maximal_faces=[(1,2,3), (1,2,4),\
          (1,2,5), (1,3,4),(1,3,5),(2,3,4),(2,3,5)])
        sage: i = 1
        sage: c = [1, 0, 0, -1, 0, 0, 1, 0, 0]
        sage: is_active(S, i, c, [2,3,6,7])
        False
        sage: c = (0, 1, 0, -1, 0, 0, 0, 0, 1)
        sage: is_active(S, i, c, [2,3,6,7])
        False

    """
    Q = simplex_laplacian_combinatorial_up(X,i)
    n = len(Q.rows())
    for j in range(n):
        if not(j in S0):
            c1 = fire_set(X, i, c, [j])
       	    b = []
	    for k in range(n):
	        if not(k in S0):
	            b = b+[c1[k]<0]
	    if True in b:
	        return False
    return True
    
def fire_set(X, i, c, S0):
    """
    Returns the configuration after firing the active vertex v along S0.

    INPUT:
     X - a simplicial complex
     c - a configuration (implemented as a list having the same length 
          as the number of cols of the i-th boundary matrix of X)
     S0 - a set of $i$-faces of X

    EXAMPLES:
        sage: i = 1
        sage: c = [1,2,3,4,5,6,7,8,9]
        sage: S0 = [2,3,6,7]
        sage: X = SimplicialComplex(maximal_faces=[(1,2,3), (1,2,4),\
            (1,2,5), (1,3,4),(1,3,5),(2,3,4),(2,3,5)])
        sage: fire_set(X, i, c, S0)
        (2, 3, 2, 3, 5, 8, 6, 7, 11)
        sage: c = [1, 0, 0, -1, 0, 0, 1, 0, 0]
        sage: fire_set(S, i, c, [2,5,7])
        (1, 0, 0, -1, 0, 0, 1, 0, 0)
        sage: fire_set(S, i, c, [3,6,8])
        (1, 0, 0, -1, 0, 0, 1, 0, 0)
        sage: c = [2, 0, 0, 4, 1, 1, 4, 2, 0] # stable wrt [2,3,6,7]
        sage: fire_set(S, i, c, [3])
        (3, 1, 0, 2, 1, 2, 4, 2, 1)
        sage: fire_set(S, i, c, [3,0,])
        (1, 0, 0, 3, 1, 3, 5, 2, 1)
        sage: fire_set(S, i, c, [3,0])
        (1, 0, 0, 3, 1, 3, 5, 2, 1)
        sage: fire_set(S, i, c, [3,0,5]) # thus c is recurrent
        (2, 0, 0, 4, 1, 1, 4, 2, 0)

    """
    Q = simplex_laplacian_combinatorial_up(X,i)
    #n = len(Q.rows())
    c1 = copy(c)
    vc = vector(c1)
    for j in S0:
        #print j, vc, Q.rows()[j]
        vc = vc - Q.rows()[j]
    return vc

def small_active_configurations(X, i, S0):
    """
    Returns True if the configuration c is ``small'' (in [0,..,d-1]^n, where
    d is the largest component of the f-vector of X) and active (ie, each 
    non-source can legally fire), False otherwise.

     INPUT:
     X - a simplicial complex
     S0 - an i-diml spanning tree of X, as a list of i-faces of X (the "sources")


    EXAMPLES:
        sage: S = SimplicialComplex(maximal_faces=[(1,2,3), (1,2,4),\
          (1,2,5), (1,3,4),(1,3,5),(2,3,4),(2,3,5)])
        sage: i = 1
        sage: L = small_active_configurations(S, i, [2,3,6,7]) # 3 days
        sage: len(L)
        1908125
        sage: 1908125.0/1953125 # prob a small config is active
        0.976960000000000

    """
    Q = simplex_laplacian_combinatorial_up(X,i)
    n = len(Q.rows())
    d = S.f_vector()[i] # bound on degree
    R = IntegerModRing(d)^n
    L = []
    for c in R:
        c = [ZZ(a) for a in c]
	if is_stable(X, i, c, S0):
	    L.append(c)
    return L
	
def is_stable(X, i, c, S0):
    """
    Returns True if the configuration c is stable (ie, no
    non-source can legally fire), False otherwise.
    
    INPUT:
     X - a simplicial complex
     c - a configuration (implemented as a list having the same length 
          as the number of cols of the i-th boundary matrix of X)
     S0 - a set of source $i$-faces of X (ie, a spanning tree)


    EXAMPLES:
        sage: S = SimplicialComplex(maximal_faces=[(1,2,3), (1,2,4),\
          (1,2,5), (1,3,4),(1,3,5),(2,3,4),(2,3,5)])
        sage: i = 1
        sage: c = [1, 0, 0, -1, 0, 0, 1, 0, 0]
        sage: is_stable(S, i, c, [2,3,6,7])
        True
        sage: c = (0, 1, 0, -1, 0, 0, 0, 0, 1)
        sage: is_stable(S, i, c, [2,3,6,7])
        True
        sage: is_stable(S, i, [1]*9, S0)
        True
        sage: is_stable(S, i, [2]*9, S0)
        False

    """
    Q = simplex_laplacian_combinatorial_up(X,i)
    n = len(Q.rows())
    fire = []
    for j in range(n):
        if not(j in S0):
            c1 = fire_set(X, i, c, [j])
       	    b = []
	    for k in range(n):
	        if not(k in S0):
	            b.append(c1[k]<0)
                    #print k, c, b, c1   
	    fire.append(True in b)
	    #print j, "fire = ", fire
    if False in fire:
        return False
    else:
        return True

def stable_configurations(X, i, S0):
    """
    Returns the set of stable configurations.

     INPUT:
     X - a simplicial complex
     S0 - a set of source $i$-faces of X (ie, a spanning tree)


    EXAMPLES:
        sage: S = SimplicialComplex(maximal_faces=[(1,2,3), (1,2,4),\
          (1,2,5), (1,3,4),(1,3,5),(2,3,4),(2,3,5)])
        sage: i = 1
        sage: L = stable_configurations(S, i, [2,3,6,7]) # very long time
        sage: l0 = len(L); l0
        268125
        sage: 268125.0/5^9 # prob a "small" config is stable
        0.137280000000000
        sage: j0 = int(random()*l0); j0
        98122
        sage: c = L[j0]; c
        [2, 0, 0, 4, 1, 1, 4, 2, 0]
        sage: fire_set(S, i, c, [1])
        (1, -3, 1, 5, 1, 1, 3, 1, 1)

    """
    Q = simplex_laplacian_combinatorial_up(X,i)
    n = len(Q.rows())
    d = S.f_vector()[i] # bound on degree
    R = IntegerModRing(d)^n
    L = []
    for c in R:
        c = [ZZ(a) for a in c]
	if is_stable(X, i, c, S0):
	    L.append(c)
	    #print c, len(L)
    return L
	
def is_recurrent(X, i, c, F, S0):
    """
    Returns True if the configuration c is returned by S0-legally firing the
    i-faces in F, False otherwise.
    
    INPUT:
     X - a simplicial complex
     c - a configuration (implemented as a list having the same length 
          as the number of cols of the i-th boundary matrix of X)
     S0, F - a set of source $i$-faces of X 


    EXAMPLES:
        sage: S = SimplicialComplex(maximal_faces=[(1,2,3), (1,2,4),\
          (1,2,5), (1,3,4),(1,3,5),(2,3,4),(2,3,5)])
        sage: i = 1
        sage: c = [1, 0, 0, -1, 0, 0, 1, 0, 0]
        sage: is_recurrent(S, i, c, [3, 0, 5], [2,3,6,7])
        True

    The number of S0-stable configurations of S is $268125 = 3 \cdot 5^4 \cdot 11\cdot 13$, S0 = [2,3,6,7].
    The number of those which are also S0-recurrent configurations for the firing set [3, 0, 5] 
    is $146250 = 2 \cdot 3^2 \cdot 5^4 \cdot 13$.

    """
    Q = simplex_laplacian_combinatorial_up(X,i)
    n = len(Q.rows())
    rec = []
    c1 = c
    for j in F:
            c1 = fire_set(X, i, c1, [j])
       	    b = []
	    for k in range(n):
	        if not(k in S0):
	            b.append(c1[k]<0)
                    #print k, c, b, c1   
	    rec.append(True in b)
	    #print j, "rec = ", rec
    if True in rec:
        return False
    else:
        return True

def is_recurrent2(X, i, c, S0):
    """
    Returns True if the configuration c+d is returned by S0-legally 
    firing all the i-faces in X, False otherwise.
    
    INPUT:
     X - a simplicial complex
     c - a configuration (implemented as a list having the same length 
          as the number of cols of the i-th boundary matrix of X)
     S0 - a set of source $i$-faces of X 


    EXAMPLES:
        sage: S = SimplicialComplex(maximal_faces=[(1,2,3), (1,2,4),\
          (1,2,5), (1,3,4),(1,3,5),(2,3,4),(2,3,5)])
        sage: i = 1
        sage: c = [1, 0, 0, -1, 0, 0, 1, 0, 0]
        sage: is_recurrent2(S, i, c, [2,3,6,7])
        True

    """
    Q = simplex_laplacian_combinatorial_up(X,i)
    n = len(Q.rows())
    rec = []
    c1 = c
    F = [j for j in range(n) if not(j in S0)]
    for j in S0:
            c1 = fire_set(X, i, c1, [j])
       	    b = []
	    for k in range(n):
	        if not(k in S0):
	            b.append(c1[k]<0)
                    #print k, c, b, c1   
	    rec.append(True in b)
	    #print j, "rec = ", rec
    for j in F:
            c1 = fire_set(X, i, c1, [j])
       	    b = []
	    for k in range(n):
	        if not(k in S0):
	            b.append(c1[k]<0)
                    #print k, c, b, c1   
	    rec.append(True in b)
	    #print j, "rec = ", rec	    
    if True in rec:
        return False
    else:
        return True