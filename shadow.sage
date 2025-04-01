def slope_to_modp(slopelist,p):
    modp = []
    for x in slopelist:
        try:
            modp.append(GF(p)(x))
        except ValueError:
            modp.append(oo)
    return modp

def check_p_divisibility(E,d,Plistgen,p):
    """
    Checks if Plistgen is p-divisible in E^D(Qp)
    """
    ED = E.quadratic_twist(d)
    EDs = ED.short_weierstrass_model()
    R = Plistgen
    R = EDs(R)
    if R[0].denominator().sqrt().valuation(p) == 0:
        try:
            EDsF = EDs.change_ring(GF(p))
            RF = EDsF(R); RF
            n = RF.additive_order()
        except ArithmeticError:
            EDF = ED.change_ring(GF(p))
            f_from_EDs_to_ED = EDs.isomorphism_to(ED)
            R_on_ED = f_from_EDs_to_ED(R)
            n = EDF(R_on_ED).additive_order()
        if n%p == 0:
            return False
    else:
        n = 1
    try:
        EDsF = EDs.change_ring(GF(p))
        K = Qp(p,100)
        EDsK = EDs.change_ring(K)
        Q = n*R

        a = K(-Q[0]/Q[1])
        F = EDsK.formal_group()
        if (F.log()(a)).valuation() > 1:
            print('R is p-divisible!')
            E._has_pdivisible_R = True
            return True
        else:
            EDsFratpts = EDsF.rational_points()
            for x in EDsFratpts:
                if x.additive_order() == p:
                    Z = x
                    break
            try:
                Aprime = EDsK.lift_x(ZZ(Z[0]))
            except UnboundLocalError:
                return False
            pAprime = p*Aprime
            if (F.log()(-pAprime[0]/pAprime[1])).valuation() == 1:
                print('R is p-divisible!')
                E._has_pdivisible_R = True
                return True
            else:
                return False
    except ArithmeticError:
        K = Qp(p,100)
        EDK = ED.change_ring(K)
        f_from_EDs_to_ED = EDs.isomorphism_to(ED)
        R_on_ED = f_from_EDs_to_ED(R)
        Q_on_ED = n*R_on_ED
        a = K(-Q_on_ED[0]/Q_on_ED[1])
        F = EDK.formal_group()
        if (F.log()(a)).valuation() > 1:
            print('R is p-divisible!')
            E._has_pdivisible_R = True
            return True
        else:
            EDFratpts = EDF.rational_points()
            for x in EDFratpts:
                if x.additive_order() == p:
                    Z = x
                    break
            try:
                Aprime = EDK.lift_x(ZZ(Z[0]))
                pAprime = p*Aprime
                if (F.log()(-pAprime[0]/pAprime[1])).valuation() == 1:
                    print('R is p-divisible!')
                    E._has_pdivisible_R = True
                    return True
                else:
                    return False
            except UnboundLocalError:
                return False


def model_and_gens_rank3_over_K(E, d, Plistgen,p):
    """
    E: elliptic curve with analytic rank 2 over Q
    d: quadratic imaginary discriminant such that E(K) has rank 3, where K = Q(sqrt(d))
    returns a minimal model of E and the 3 generators of E(K)
    """
    assert E.rank() == 2
    assert E.analytic_rank() == 2
    Ed = E.quadratic_twist(d)
    Eds = Ed.short_weierstrass_model()
    Eds_gen = Eds(Plistgen)
    Qp200 = Qp(p,200)
    EdQp = Eds.change_ring(Qp200)
    EdgQp = EdQp(Eds_gen[0],Eds_gen[1])
    check_p_divisibility(E,d,Plistgen,p)
    Enew = EllipticCurve([Eds.a4()/d^2, Eds.a6()/d^3])
    K.<D> = QuadraticField(d)
    h = K.class_number()
    print('h_K:', h)
    if h%p == 0:
        print('p divides h_K !')
    EnewK = Enew.change_ring(K)
    EnewKpoints = [EnewK(P) for P in Enew.gens()] + [EnewK(Eds_gen[0]/D^2, (Eds_gen[1]/D^3))]
    EnewKmin = EnewK.change_ring(QQ).minimal_model().change_ring(K)
    f = EnewK.isomorphism_to(EnewKmin)
    EnewKminpoints = [f(P) for P in EnewKpoints]
    return EnewKmin, EnewKminpoints

def embeddings(K,p,prec):
    """
    The embedding(s) $K=\Q(\sqrt(D)) \into \Q_p$.
    """
    Q = Qp(p,prec)
    OK = K.maximal_order()
    pOK = factor(p*OK)
    if (len(pOK) == 2 and pOK[0][1] == 1):
        R = Q['x']
        r1, r2 = R(K.defining_polynomial()).roots()
        psi1 = K.hom([r1[0]])
        psi2 = K.hom([r2[0]])
        return [psi1, psi2]
    else:
        F = Q.extension(K.defining_polynomial(),names='a')
        a = F.gen()
        psi = self._psis = [K.hom([a])]
        return psi

def anticyc_padic_height(E, P, d, p, prec, multi, timelimit):
    """
    computes the anticyclotomic p-adic height of P
    """
    Qp = pAdicField(p,prec)
    EQp = E.change_ring(Qp)
    assert E.is_good(p) and E.is_ordinary(p)
    K.<D> = QuadraticField(d)
    h = K.class_number()
    if h%p == 0:
        E._has_pdivisible_hK = True
    OK = K.maximal_order()

    assert len(factor(p*OK)) == 2

    pi, pibar = factor(p*OK)
    pi = pi[0]
    pibar = pibar[0]

    p1 = (pi^h).gens_reduced()[0]
    p2 = (pibar^h).gens_reduced()[0]

    psi1, psi2 = embeddings(K,p,prec)
    #embedding is chosen so that p1 is 'pi'
    if psi1(p1).valuation() > 0:
        embedding = psi1
    else:
        embedding = psi2
    assert embedding(p1).valuation() > 0
    m0 = lcm(E.tamagawa_numbers())
    R = multi*m0*P
    try_orders = divisors(len(E.change_ring(GF(p)).rational_points()))
    for ord in try_orders:
        B1 = ord*R
        B1conj = (B1[0].conjugate(), B1[1].conjugate())
        B1 = EQp(embedding(B1[0]), embedding(B1[1]))
        B1conj = EQp(embedding(B1conj[0]), embedding(B1conj[1]))
        tB1 = -B1[0]/B1[1]
        tB1conj = -B1conj[0]/B1conj[1]
        if tB1.valuation() > 0 and tB1conj.valuation() > 0 and B1[0].valuation() < 0 and B1[1].valuation() <0 and B1conj[0].valuation() < 0 and B1conj[1].valuation() < 0:
            n = ord
            break
    T = n*R

    if n%2 == 1:
        fn = E.change_ring(QQ).division_polynomial(n)
    else:
        fn = E.change_ring(QQ).division_polynomial(n, two_torsion_multiplicity = 1)

    xR = R[0]
    xT = T[0]
    if n%2 == 1:
        fn_at_R = fn(xR)
    else:
        fn_at_R = fn((R[0], R[1]))

    if xR == 0:
        d_at_R = 1
    else:
        d_at_R = K.ideal(xR).denominator()
    alarm(timelimit)

    try:
        d_at_R = prod([a^(ZZ(e/2)) for (a,e) in factor(d_at_R)])
        d_at_R = (d_at_R^h).gens_reduced()[0]
    except AttributeError:
        assert (h == 1 or d_at_R == K.ideal(1))
        d_at_R = prod([a^(ZZ(e/2)) for (a,e) in factor(d_at_R)])
    alarm(timelimit)
    d_of_nR_to_h = d_at_R^(n^2)*fn_at_R^h
    dbydconjnR = d_of_nR_to_h/d_of_nR_to_h.conjugate()
    value_away_from_p_nR_for_embedding = K.ideal((d_of_nR_to_h.conjugate()/d_of_nR_to_h)).gens_reduced()[0]
    height_away_from_p_nR_via_embedding = 1/(p*h)*log(embedding(value_away_from_p_nR_for_embedding),0)  #this is the height away from p of T=nR
    sig = E.change_ring(QQ).padic_sigma(p,prec)
    Tconj = (T[0].conjugate(), T[1].conjugate())
    height_at_p_nR  = 1/(p)*log(sig(embedding(-T[0]/T[1]))/sig(embedding(-Tconj[0]/Tconj[1])),0) #this is the height at p of T = nR
    E._at_p = height_at_p_nR
    height_T = height_at_p_nR + height_away_from_p_nR_via_embedding
    height_R = 1/n^2*(height_T)
    height_P = 1/(m0^2*multi^2)*height_R
    return height_P

def shadow_line_slope(E, d, p, frommagma, prec, timelimit, mwgens):
    #slopes are truncated at O(p^50)
    Enew,pts = model_and_gens_rank3_over_K(E, d, frommagma,p)
    pt1, pt2, Q = pts
    if mwgens != None:
        pt1 = mwgens[0]
        pt1 = Enew(pt1[0],pt1[1])
        pt2 = mwgens[1]
        pt2 = Enew(pt2[0],pt2[1])
    print("pt1: ", pt1)
    print("pt2: ", pt2)
    Ep = E.change_ring(GF(p))
    P1 = E(pt1[0],pt1[1])
    P2 = E(pt2[0],pt2[1])
    n = E.Np(p)
    P1p = Ep(P1)
    P2p = Ep(P2)
    b1 = P1p.additive_order().valuation(p)
    b2 = P2p.additive_order().valuation(p)
    if b1 == 1 and b2 == 1:
        for c in range(1,p):
            if (P1p +c*P2p).additive_order().valuation(p) == 0:
                break
        case = 'case3'
        #then generators of H are P1 + c*P2, p*P2
    elif b1 == 0 and b2 == 0:
        case = 'case1'
        #then generators of H are same as before
    else:
        case = 'case2'
        cval = [p^b1,p^b2]
        #then generators of H are p^b1*P1, p^b2*P2
    if case == 'case3':
        G1 = pt1 + c*pt2 + Q
        G2 = pt2 + Q #need to multiply height by 1/p, since this should be p*pt2
    else:
        G1 = pt1 + Q
        G2 = pt2 + Q
    sys.stdout.flush()
    hpt1 = anticyc_padic_height(E, G1, d, p, prec, 1, timelimit)
    hpt2 = anticyc_padic_height(E, G2, d, p, prec, 1, timelimit)
    if case == 'case1':
        C =  -hpt1/hpt2 +O(p^50)  ##this is (1,p) *1/p or (p,1) *p for case 2
    elif case == 'case3':
        C = -hpt1/(p*hpt2) +O(p^50)
    else:
        C = -cval[0]*hpt1/(cval[1]*hpt2) + O(p^50)
    return C


def partp(modlist, p):
    listminusp = [modlist.count(i) for i in range(p)]
    listinfty = len(modlist) - sum(listminusp)
    return listminusp + [listinfty]

def shadow_distribution(E, p, Dlist, DPlist, timelimit, mwgens = None):
    print('Elliptic curve E: %s (Cremona label)'%E)
    sys.stdout.flush()
    try:
        E = EllipticCurve(E)
    except:
        pass
    print("List of discriminants: ", Dlist)
    print(50*"*")
    slopelist = []
    sys.stdout.flush()
    for i in range(len(Dlist)):
        if Dlist[i] < 0:
            try:
                if E.Np(p)%p == 0:
                    sys.stdout.flush()
                print('example#: ', i+1)
                print("D = ", Dlist[i])
                sys.stdout.flush()
                A = shadow_line_slope(E,Dlist[i],p,DPlist[i], 500, timelimit, mwgens)
                Lslope = A%p
                print('slope: ', A)
                sys.stdout.flush()
                slopelist.append(A)
                print('stats: ', partp(slope_to_modp(slopelist,p),p))
                sys.stdout.flush()
                print(20*'-')
            except RuntimeError:
                print('RuntimeError at %s'%(Dlist[i]))
                print(20*'-')
                sys.stdout.flush()
                pass
            except AssertionError:
                print('AssertionError at %s'%(Dlist[i]))
                print(20*'-')
                sys.stdout.flush()
                pass
            except (TypeError,ValueError):
                print('ValueError at %s'%(Dlist[i]))
                print('retrying with increased precision')
                A = shadow_line_slope(E,Dlist[i],p,DPlist[i], 1300, timelimit, mwgens)
                print('slope: ', A)
                slopelist.append(A)
                print('stats: ', partp(slope_to_modp(slopelist,p),p))
                print(20*'-')
                sys.stdout.flush()
            except AlarmInterrupt:
                print('AlarmInterrupt at %s'%(Dlist[i]))
                print(20*'-')
                sys.stdout.flush()
                pass

    print('slope list: ', slopelist)
    sys.stdout.flush()
