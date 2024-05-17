import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import networkx as nx


# Define your build_trinomial_tree_vasicek function
def build_trinomial_tree_vasicek(delta_t, T_length, r_0, alpha, beta, K, f, k_v, theta_v, sigma_v):
    # Initialize the tree with the r_0 node
    t = np.arange(0, T_length + delta_t, delta_t)
    T = np.arange(0, T_length + 1)
    V = np.sqrt(sigma_v**2/(2*theta_v)*(1-np.exp(-2*theta_v*delta_t)))
    delta=np.sqrt(3)*V
    x = [[r_0]] 
    # Generate subsequent T[-1],
    for i in range(1, len(t)):
        # Number of nodes in the current layer: 1, 3, 5, 7, ...
        num_nodes = 2 * i + 1
        # Create the current layer
        current_layer = []
        # Calculate the starting position to center the nodes
        start_pos = -(num_nodes - 1) * delta / 2 + r_0
        for j in range(num_nodes):
            # Calculate the position of each node
            pos = start_pos + j * delta
            current_layer.append(pos)
        # Add the current layer to the x
        x.append(current_layer)
    M = []
    for i in range(0, len(t)):
        M_i = []
        for j in range(len(x[i])):
            M_ij = x[i][j]*np.exp(-theta_v*0.5)+k_v/theta_v*(1-np.exp(-theta_v*delta_t))
            M_i.append(M_ij)
        M.append(M_i)
    print('M',M)
    print('x',x)
    print('V',V)
    # print('p_u',p_u)
    # print('p_u',p_u)
    # print('p_u',p_u)r=

    p_u=[]
    p_m=[]
    p_d=[]
    for i in range(0, len(t)-1):
        num_nodes = 2 * i + 1
        current_p_u = []
        current_p_m = []
        current_p_d = []
        start_pos = -(num_nodes - 1) * delta / 2 + r_0
        for j in range(num_nodes):
            k=j+1
            p_u_pos =1/6+(M[i][j]-x[i+1][k])**2/(6*V**2)+(M[i][j]-x[i+1][k])/(2*np.sqrt(3)*V)
            p_m_pos =2/3-(M[i][j]-x[i+1][k])**2/(3*V**2)
            p_d_pos =1/6+(M[i][j]-x[i+1][k])**2/(6*V**2)-(M[i][j]-x[i+1][k])/(2*np.sqrt(3)*V)
            current_p_u.append(p_u_pos)
            current_p_m.append(p_m_pos)
            current_p_d.append(p_d_pos)
        p_u.append(current_p_u)
        p_m.append(current_p_m)
        p_d.append(current_p_d)
    print('p_u',p_u)
    print('p_m',p_m)
    print('p_d',p_d)
    
    # P_market=np.exp(-np.array(f)/100)
    cumsum_f = np.cumsum(f)
    cumsum_f = np.insert(cumsum_f, 0, 0)
    P_market = np.exp(-cumsum_f*delta_t / 100)
    print("P_market", P_market)
    phi=[]
    Q=[[1]]
    for i in range(0, len(t)):
        summa = 0
        for h in range(len(Q[i])):
            Q_ih = Q[i][h]
            summa += Q_ih * np.exp(-(delta_t)*x[i][h])
        # print(P_market)
        # print(P_market[i+1])
        # print(t[i+1])
        # print("summa", summa)
        # print("P_market[i+1]", P_market[i+1])
        phi_i = (np.log(summa)-np.log(P_market[i+1]))/(delta_t)
        phi.append(phi_i)
        if i < len(t)-1:
            Q_plus = np.zeros(len(x[i+1]))
            for j in range(len(x[i])):
                Q_plus[j] +=  Q[i][j]*p_d[i][j]*np.exp(-(t[i+1]-t[i])*(x[i][j]+phi[i]))
                Q_plus[j+1] +=  Q[i][j]*p_m[i][j]*np.exp(-(t[i+1]-t[i])*(x[i][j]+phi[i]))       
                Q_plus[j+2] +=  Q[i][j]*p_u[i][j]*np.exp(-(t[i+1]-t[i])*(x[i][j]+phi[i]))
            Q.append(list(Q_plus))
            #  Q_plush =  Q[i][h]*p[i][h]
    print("Q",Q)
    print("Phi",phi)
    print("x",x)
    for i in range(0, len(t)):
        for j in range(0, len(x[i])):
            x[i][j] += phi[i]
    P = []
    r=x
    print("r", r)
    for Ts in range(0,len(T)):
        PTs=[]
        expiry_time = np.where(t == Ts)[0][0]
        for i in range(expiry_time, -1, -1):
            # print(expiry_time)

            if expiry_time == i:
                PTs.append(list(np.ones(len(x[i]))))
            else:
                current_P=[]
                for j in range(len(x[i])):
                    # print(Ts,  expiry_time,len(x[i]),i, j,PTs)

                    P_pos = (p_u[i][j]*PTs[-1][j+2]+p_m[i][j]*PTs[-1][j+1]+p_d[i][j]*PTs[-1][j])*np.exp(-r[i][j]*(delta_t))
                    current_P.append(P_pos)
                PTs.append(current_P)
        PTs.reverse()
        P.append(PTs)
    tau=np.diff(T)
    # print(tau)
    k=len(T)-1
    IRS=[]
    t_alpha=np.where(t == alpha)[0][0]
    t_beta=np.where(t == beta)[0][0]
    
    for l in range(alpha, beta+1):
        t_l=np.where(t == l)[0][0]
        IRS_l=[]
        for j in range(len(x[t_l])):
            # print(l,i,len(x[l]),len(T),P[k])
            sum=0
            for s in range(l+1,k+1):
                # print(sum,l,k,s,j)
                # print(P[s][t_l])
                sum=sum+tau[s-1]*K*P[s][t_l][j]
                
            IRS_lj=1-P[k][t_l][j]-sum
            IRS_l.append(IRS_lj)
        IRS.append(IRS_l)
    IRS_O=[]
    for l in range(beta,-1,-1):
        if l== beta:
            print(IRS)
            IRS_O.append(np.maximum(IRS[-1],np.zeros(len(IRS[-1]))))
        elif l>=alpha:
            # for j in range(len(x[l])):
            # print(len(x[l]), l,IRS_O[-1][j+2],p_u[l][j])
            T_l1=np.where(t == T[l+1])[0][0]
            T_l=np.where(t == T[l])[0][0]
            EV=[]
            # print(T_l1-T_l)
            for n in range(T_l1,T_l-1,-1):
                if n==T_l1:
                    EV.append(IRS_O[-1])
                
                else:
                    EV_n=[]
                    for m in range(len(x[n])):
                        # print(p_u[n][m])
                        # print(IRS_O)
                        # print(n,m,l,len(x[n]))
                        # print(EV)
                        # print(EV[-1][m+2])
                        EV_nm=(p_u[n][m]*EV[-1][m+2]+p_m[n][m]*EV[-1][m+1]+p_d[n][m]*EV[-1][m])*np.exp(-(t[n+1]-t[n])*r[n][m])
                        EV_n.append(EV_nm)
                    EV.append(EV_n)
            O_P=np.maximum(EV[-1],IRS[T[l]-alpha])
            IRS_O.append(O_P)
        else:
            T_l1=np.where(t == T[l+1])[0][0]
            T_l=np.where(t == T[l])[0][0]
            EV=[]
            # print(T_l1-T_l)
            for n in range(T_l1,T_l-1,-1):
                if n==T_l1:
                    EV.append(IRS_O[-1])
                
                else:
                    EV_n=[]
                    for m in range(len(x[n])):
                        # print(p_u[n][m])
                        # print(IRS_O)
                        # print(n,m,l,len(x[n]))
                        # print(EV)
                        # print(EV[-1][m+2])
                        EV_nm=(p_u[n][m]*EV[-1][m+2]+p_m[n][m]*EV[-1][m+1]+p_d[n][m]*EV[-1][m])*np.exp(-(t[n+1]-t[n])*r[n][m])
                        EV_n.append(EV_nm)
                    EV.append(EV_n)
            O_P=np.array(EV[-1])
            IRS_O.append(O_P)
    price=IRS_O[-1][0]
                # print(EV_n)
            # EV = EV_n[-1] r = 
            
            # EV_lj=p_u[l][j]*IRS_O[-1][j+2]+p_m[l][j]*IRS_O[-1][j+1]+p_d[l][j]*IRS_O[-1][j]
                
    return r, x, p_u, p_m, p_d, P, IRS, IRS_O, price

def build_trinomial_tree_cir(delta_t, T_length, r_0, alpha, beta, K, f, alpha_cir, beta_cir, sigma_cir):
    # Initialize the tree with the r_0 node
    t = np.arange(0, T_length + delta_t, delta_t)
    T = np.arange(0, T_length + 1)
    V = sigma_cir/2
    delta=np.sqrt(3)*V
    x = [[r_0]] 
    # Generate subsequent T[-1],
    for i in range(1, len(t)):
        # Number of nodes in the current layer: 1, 3, 5, 7, ...
        num_nodes = 2 * i + 1
        # Create the current layer
        current_layer = []
        # Calculate the starting position to center the nodes
        start_pos = -(num_nodes - 1) * delta / 2 + r_0
        for j in range(num_nodes):
            # Calculate the position of each node
            pos = start_pos + j * delta
            current_layer.append(pos)
        # Add the current layer to the x
        x.append(current_layer)
    M = []
    for i in range(0, len(t)):
        M_i = []
        for j in range(len(x[i])):
            M_ij = x[i][j]+((alpha_cir/2-sigma_cir**2/8)/x[i][j]-beta_cir/2*x[i][j])*delta_t
            M_i.append(M_ij)
        M.append(M_i)
    print('M',M)
    print('x',x)
    print('V',V)
    # print('p_u',p_u)
    # print('p_u',p_u)
    # print('p_u',p_u)r=

    p_u=[]
    p_m=[]
    p_d=[]
    for i in range(0, len(t)-1):
        num_nodes = 2 * i + 1
        current_p_u = []
        current_p_m = []
        current_p_d = []
        start_pos = -(num_nodes - 1) * delta / 2 + r_0
        for j in range(num_nodes):
            k=j+1
            p_u_pos =1/6+(M[i][j]-x[i+1][k])**2/(6*V**2)+(M[i][j]-x[i+1][k])/(2*np.sqrt(3)*V)
            p_m_pos =2/3-(M[i][j]-x[i+1][k])**2/(3*V**2)
            p_d_pos =1/6+(M[i][j]-x[i+1][k])**2/(6*V**2)-(M[i][j]-x[i+1][k])/(2*np.sqrt(3)*V)
            current_p_u.append(p_u_pos)
            current_p_m.append(p_m_pos)
            current_p_d.append(p_d_pos)
        p_u.append(current_p_u)
        p_m.append(current_p_m)
        p_d.append(current_p_d)
    print('p_u',p_u)
    print('p_m',p_m)
    print('p_d',p_d)

    x=[[element ** 2 for element in sublist] for sublist in x]

    cumsum_f = np.cumsum(f)
    cumsum_f = np.insert(cumsum_f, 0, 0)
    P_market = np.exp(-cumsum_f*delta_t / 100)
    print("P_market", P_market)
    phi=[]
    Q=[[1]]
    for i in range(0, len(t)):
        summa = 0
        for h in range(len(Q[i])):
            Q_ih = Q[i][h]
            summa += Q_ih * np.exp(-(delta_t)*x[i][h])
        # print(P_market)
        # print(P_market[i+1])
        # print(t[i+1])
        # print("summa", summa)
        # print("P_market[i+1]", P_market[i+1])
        phi_i = (np.log(summa)-np.log(P_market[i+1]))/(delta_t)
        phi.append(phi_i)
        if i < len(t)-1:
            Q_plus = np.zeros(len(x[i+1]))
            for j in range(len(x[i])):
                Q_plus[j] +=  Q[i][j]*p_d[i][j]*np.exp(-(t[i+1]-t[i])*(x[i][j]+phi[i]))
                Q_plus[j+1] +=  Q[i][j]*p_m[i][j]*np.exp(-(t[i+1]-t[i])*(x[i][j]+phi[i]))       
                Q_plus[j+2] +=  Q[i][j]*p_u[i][j]*np.exp(-(t[i+1]-t[i])*(x[i][j]+phi[i]))
            Q.append(list(Q_plus))
            #  Q_plush =  Q[i][h]*p[i][h]
    print("Q",Q)
    print("Phi",phi)
    print("x",x)
    for i in range(0, len(t)):
        for j in range(0, len(x[i])):
            x[i][j] += phi[i]
    P = []
    r=x
    print("r", r)
    for Ts in range(0,len(T)):
        PTs=[]
        expiry_time = np.where(t == Ts)[0][0]
        for i in range(expiry_time, -1, -1):
            # print(expiry_time)

            if expiry_time == i:
                PTs.append(list(np.ones(len(x[i]))))
            else:
                current_P=[]
                for j in range(len(x[i])):
                    # print(Ts,  expiry_time,len(x[i]),i, j,PTs)

                    P_pos = (p_u[i][j]*PTs[-1][j+2]+p_m[i][j]*PTs[-1][j+1]+p_d[i][j]*PTs[-1][j])*np.exp(-r[i][j]*(delta_t))
                    current_P.append(P_pos)
                PTs.append(current_P)
        PTs.reverse()
        P.append(PTs)
    tau=np.diff(T)
    # print(tau)
    k=len(T)-1
    IRS=[]
    t_alpha=np.where(t == alpha)[0][0]
    t_beta=np.where(t == beta)[0][0]
    
    for l in range(alpha, beta+1):
        t_l=np.where(t == l)[0][0]
        IRS_l=[]
        for j in range(len(x[t_l])):
            # print(l,i,len(x[l]),len(T),P[k])
            sum=0
            for s in range(l+1,k+1):
                # print(sum,l,k,s,j)
                # print(P[s][t_l])
                sum=sum+tau[s-1]*K*P[s][t_l][j]
                
            IRS_lj=1-P[k][t_l][j]-sum
            IRS_l.append(IRS_lj)
        IRS.append(IRS_l)
    IRS_O=[]
    for l in range(beta,-1,-1):
        if l== beta:
            print(IRS)
            IRS_O.append(np.maximum(IRS[-1],np.zeros(len(IRS[-1]))))
        elif l>=alpha:
            # for j in range(len(x[l])):
            # print(len(x[l]), l,IRS_O[-1][j+2],p_u[l][j])
            T_l1=np.where(t == T[l+1])[0][0]
            T_l=np.where(t == T[l])[0][0]
            EV=[]
            # print(T_l1-T_l)
            for n in range(T_l1,T_l-1,-1):
                if n==T_l1:
                    EV.append(IRS_O[-1])
                
                else:
                    EV_n=[]
                    for m in range(len(x[n])):
                        # print(p_u[n][m])
                        # print(IRS_O)
                        # print(n,m,l,len(x[n]))
                        # print(EV)
                        # print(EV[-1][m+2])
                        EV_nm=(p_u[n][m]*EV[-1][m+2]+p_m[n][m]*EV[-1][m+1]+p_d[n][m]*EV[-1][m])*np.exp(-(t[n+1]-t[n])*r[n][m])
                        EV_n.append(EV_nm)
                    EV.append(EV_n)
            O_P=np.maximum(EV[-1],IRS[T[l]-alpha])
            IRS_O.append(O_P)
        else:
            T_l1=np.where(t == T[l+1])[0][0]
            T_l=np.where(t == T[l])[0][0]
            EV=[]
            # print(T_l1-T_l)
            for n in range(T_l1,T_l-1,-1):
                if n==T_l1:
                    EV.append(IRS_O[-1])
                
                else:
                    EV_n=[]
                    for m in range(len(x[n])):
                        # print(p_u[n][m])
                        # print(IRS_O)
                        # print(n,m,l,len(x[n]))
                        # print(EV)
                        # print(EV[-1][m+2])
                        EV_nm=(p_u[n][m]*EV[-1][m+2]+p_m[n][m]*EV[-1][m+1]+p_d[n][m]*EV[-1][m])*np.exp(-(t[n+1]-t[n])*r[n][m])
                        EV_n.append(EV_nm)
                    EV.append(EV_n)
            O_P=np.array(EV[-1])
            IRS_O.append(O_P)
    price=IRS_O[-1][0]
                # print(EV_n)
            # EV = EV_n[-1] r = 
            
            # EV_lj=p_u[l][j]*IRS_O[-1][j+2]+p_m[l][j]*IRS_O[-1][j+1]+p_d[l][j]*IRS_O[-1][j]
                
    return r, x, p_u, p_m, p_d, P, IRS, IRS_O, price


# Streamlit app
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Define the plot_trinomial_tree function
def plot_trinomial_tree(r):
    fig, ax = plt.subplots(figsize=(14, 8))
    G = nx.DiGraph()

    pos = {}
    labels = {}
    node_count = 0

    # Map each node to its unique ID for easy edge creation
    node_map = {}
    current_id = 0
    for level, nodes in enumerate(r):
        node_map[level] = []
        for i, rate in enumerate(nodes):
            node_map[level].append(current_id)
            G.add_node(current_id, level=level, rate=rate)
            pos[current_id] = (level, -i)
            labels[current_id] = f'{rate:.2f}'
            current_id += 1

    # Create edges
    for level in range(len(r) - 1):
        for i in range(len(r[level])):
            # Current node ID
            current_node = node_map[level][i]
            # Connect to next level's nodes 
            G.add_edge(current_node, node_map[level + 1][i])
            G.add_edge(current_node, node_map[level + 1][i + 1])
            G.add_edge(current_node, node_map[level + 1][i + 2])

    nx.draw(G, pos, labels=labels, with_labels=True, node_size=700, node_color="skyblue", font_size=10, ax=ax, arrowsize=10)
    plt.title('Trinomial Tree Visualization')
    plt.xlabel('Time Steps')
    plt.ylabel('Rates')
    plt.grid(True)
    return fig

# Streamlit app
st.markdown(
    """
    <style>
    .main .block-container {
        max-width: 1600px;
        padding-left: 5%;
        padding-right: 5%;
    }
    </style>
    """,
    unsafe_allow_html=True
)
st.title('Trinomial Tree Pricing Model')
col1, col2 = st.columns(2)

def display_in_box(static_text, variable_text):
    st.markdown(f"""
    <style>
    .box {{
        border: 2px solid #4CAF50;
        padding: 20px;
        border-radius: 10px;
        background-color: #f9f9f9;
        color: #333;
    }}
    </style>
    <div class="box">
        <p>{static_text}</p>
        <p>{variable_text}</p>
    </div>
    """, unsafe_allow_html=True)

with col1:
    st.header("Vasicek")
    # Input fields
    r_0 = st.number_input('r_0', value=0.05)
    k_v = st.number_input('k', value=0.0048, step=0.001, format="%.4f")
    theta_v = st.number_input('theta', value=0.10)
    sigma_v = st.number_input('sigma_v', value=0.012, step=0.005, format="%.3f")
    K = st.number_input('K', value=0.049, step=0.005, format="%.3f")

    alpha = st.number_input('T_alpha', value=2, step=1)
    beta = st.number_input('T_beta', value=5, step=1)
    

    T_length = st.number_input('T length', value=7, step=1)
    delta_t = st.selectbox('delta_t', [1, 1/2, 1/3, 1/4], index=1)

    k_v_vals = np.arange(max(0,k_v-10*0.0005), k_v+10*0.0005, (k_v+10*0.0005-max(0,k_v-10*0.0005))/20)
    theta_v_vals = np.arange(max(0,theta_v-10*0.01), theta_v+10*0.01, (theta_v+10*0.01-max(0,theta_v-10*0.01))/20)
    sigma_v_vals = np.arange(max(0,sigma_v-10*0.001), sigma_v+10*0.001, (sigma_v+10*0.001-max(0,sigma_v-10*0.001))/20)
    K_vals = np.arange(max(0,K-10*0.005), K+10*0.005, (K+10*0.005-max(0,K-10*0.005))/20)


    # Input field for forward rate curve
    forward_rate_input = st.text_input('Forward rate curve (comma-separated)', '5, 5.1, 5.2, 5.15, 5, 4.95, 4.9, 4.8, 4.7, 4.75, 4.8, 4.9, 4.95, 5, 5.1')

    # Parse the forward rate curve input
    f = list(map(float, forward_rate_input.split(',')))

    # Ensure the forward rates array `f` matches the length of `t`
    t = np.arange(0, T_length + delta_t, delta_t)
    if len(f) != len(t):
        st.error('The number of forward rates must match the length of time steps t.')
    else:
        r, x, p_u, p_m, p_d, P, IRS, IRS_O, price = build_trinomial_tree_vasicek(delta_t, T_length, r_0, int(alpha), int(beta), K, f, k_v, theta_v, sigma_v)
        display_in_box('The price with the original parameters:', price)
        # Generate prices 1
        prices = np.zeros(len(k_v_vals))

        for i, k_v_i in enumerate(k_v_vals):
            r, x, p_u, p_m, p_d, P, IRS, IRS_O, prices[i] = build_trinomial_tree_vasicek(
                delta_t, T_length, r_0, int(alpha), int(beta), K, f, k_v_i, theta_v, sigma_v)

        # Visualization
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(k_v_vals, prices, marker='o')
        ax1.plot(k_v_vals[10], prices[10], marker='o', color='red', markersize=10)
        ax1.set_title('Price vs k')
        ax1.set_xlabel('k')
        ax1.set_ylabel('Price')
        ax1.grid(True)
        st.pyplot(fig1)

        # Generate prices 2
        prices = np.zeros(len(theta_v_vals))

        for i, theta_v_i in enumerate(theta_v_vals):
            r, x, p_u, p_m, p_d, P, IRS, IRS_O, prices[i] = build_trinomial_tree_vasicek(
                delta_t, T_length, r_0, int(alpha), int(beta), K, f, k_v, theta_v_i, sigma_v)

        # Visualization
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(theta_v_vals, prices, marker='o')
        ax1.plot(theta_v_vals[10], prices[10], marker='o', color='red', markersize=10)
        ax1.set_title('Price vs theta')
        ax1.set_xlabel('theta')
        ax1.set_ylabel('Price')
        ax1.grid(True)
        st.pyplot(fig1)

        # Generate prices 3
        prices = np.zeros(len(sigma_v_vals))

        for i, sigma_v_i in enumerate(sigma_v_vals):
            r, x, p_u, p_m, p_d, P, IRS, IRS_O, prices[i] = build_trinomial_tree_vasicek(
                delta_t, T_length, r_0, int(alpha), int(beta), K, f, k_v, theta_v, sigma_v_i)

        # Visualization
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(sigma_v_vals, prices, marker='o')
        ax1.plot(sigma_v_vals[10], prices[10], marker='o', color='red', markersize=10)
        ax1.set_title('Price vs sigma')
        ax1.set_xlabel('sigma')
        ax1.set_ylabel('Price')
        ax1.grid(True)
        st.pyplot(fig1)

        # Generate prices 4
        prices = np.zeros(len(K_vals))

        for i, K_i in enumerate(K_vals):
            r, x, p_u, p_m, p_d, P, IRS, IRS_O, prices[i] = build_trinomial_tree_vasicek(
                delta_t, T_length, r_0, int(alpha), int(beta), K_i, f, k_v, theta_v, sigma_v)

        # Visualization
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(K_vals, prices, marker='o')
        ax1.plot(K_vals[10], prices[10], marker='o', color='red', markersize=10)
        ax1.set_title('Price vs K')
        ax1.set_xlabel('K')
        ax1.set_ylabel('Price')
        ax1.grid(True)
        st.pyplot(fig1)

        # Trinomial tree visualization
        st.subheader('Trinomial Tree Visualization')
        fig2 = plot_trinomial_tree(r)
        st.pyplot(fig2)


with col2:
    st.header("CIR")
    # Input fields
    r_0 = st.number_input('r_0 ', value=0.05)
    alpha_cir = st.number_input('alpha', value=0.0048, step=0.001, format="%.4f")
    beta_cir = st.number_input('beta', value=0.10)
    sigma_cir = st.number_input('sigma_cir', value=0.0548, step=0.001, format="%.4f")
    K = st.number_input('K ', value=0.049, step=0.005, format="%.3f")
    
    alpha = st.number_input('T_alpha ', value=2, step=1)
    beta = st.number_input('T_beta ', value=5, step=1)
    

    T_length = st.number_input('T length ', value=7, step=1)
    delta_t = st.selectbox('delta_t ', [1, 1/2, 1/3, 1/4], index=1)

    # Convert K_range slider to a range of values
    alpha_cir_vals= np.arange(max(0,alpha_cir-10*0.0005), alpha_cir+10*0.0005, (alpha_cir+10*0.0005-max(0,alpha_cir-10*0.0005))/20)
    beta_cir_vals= np.arange(max(0,beta_cir-10*0.01), beta_cir+10*0.01, (beta_cir+10*0.01-max(0,beta_cir-10*0.01))/20)
    sigma_cir_vals= np.arange(max(0,sigma_v-10*0.005), sigma_v+10*0.005, (sigma_v+10*0.005-max(0,sigma_v-10*0.005))/20)

    # Generate the time steps based on T and delta_t


    # Input field for forward rate curve
    forward_rate_input = st.text_input('Forward rate curve (comma-separated) ', '5, 5.1, 5.2, 5.15, 5, 4.95, 4.9, 4.8, 4.7, 4.75, 4.8, 4.9, 4.95, 5, 5.1')

    # Parse the forward rate curve input
    f = list(map(float, forward_rate_input.split(',')))

    # Ensure the forward rates array `f` matches the length of `t`
    t = np.arange(0, T_length + delta_t, delta_t)
    if len(f) != len(t):
        st.error('The number of forward rates must match the length of time steps t.')
    else:
        r, x, p_u, p_m, p_d, P, IRS, IRS_O, price = build_trinomial_tree_cir(delta_t, T_length, r_0, int(alpha), int(beta), K, f, k_v, theta_v, sigma_v)
        display_in_box('The price with the original parameters:', price)
        # Generate prices 1
        prices = np.zeros(len(alpha_cir_vals))

        for i, alpha_cir_i in enumerate(alpha_cir_vals):
            r, x, p_u, p_m, p_d, P, IRS, IRS_O, prices[i] = build_trinomial_tree_cir(
                delta_t, T_length, r_0, int(alpha), int(beta), K, f, alpha_cir_i, beta_cir, sigma_cir)

        # Visualization
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(alpha_cir_vals, prices, marker='o')
        ax1.plot(alpha_cir_vals[10], prices[10], marker='o', color='red', markersize=10)
        ax1.set_title('Price vs alpha')
        ax1.set_xlabel('alpha')
        ax1.set_ylabel('Price')
        ax1.grid(True)
        st.pyplot(fig1)

        # Generate prices 2
        prices = np.zeros(len(beta_cir_vals))

        for i, beta_cir_i in enumerate(beta_cir_vals):
            r, x, p_u, p_m, p_d, P, IRS, IRS_O, prices[i] = build_trinomial_tree_cir(
                delta_t, T_length, r_0, int(alpha), int(beta), K, f, alpha_cir, beta_cir_i, sigma_cir)

        # Visualization
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(beta_cir_vals, prices, marker='o')
        ax1.plot(beta_cir_vals[10], prices[10], marker='o', color='red', markersize=10)
        ax1.set_title('Price vs beta')
        ax1.set_xlabel('beta')
        ax1.set_ylabel('Price')
        ax1.grid(True)
        st.pyplot(fig1)

        # Generate prices 3
        prices = np.zeros(len(sigma_cir_vals))

        for i, sigma_cir_i in enumerate(sigma_cir_vals):
            r, x, p_u, p_m, p_d, P, IRS, IRS_O, prices[i] = build_trinomial_tree_cir(
                delta_t, T_length, r_0, int(alpha), int(beta), K, f, alpha_cir, beta_cir, sigma_cir_i)

        # Visualization
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(sigma_cir_vals, prices, marker='o')
        ax1.plot(sigma_cir_vals[10], prices[10], marker='o', color='red', markersize=10)
        ax1.set_title('Price vs sigma')
        ax1.set_xlabel('sigma')
        ax1.set_ylabel('Price')
        ax1.grid(True)
        st.pyplot(fig1)

        # Generate prices 4
        prices = np.zeros(len(K_vals))

        for i, K_i in enumerate(K_vals):
            r, x, p_u, p_m, p_d, P, IRS, IRS_O, prices[i] = build_trinomial_tree_cir(
                delta_t, T_length, r_0, int(alpha), int(beta), K_i, f, alpha_cir, beta_cir, sigma_cir)

        # Visualization
        fig1, ax1 = plt.subplots(figsize=(10, 6))
        ax1.plot(K_vals, prices, marker='o')
        ax1.plot(K_vals[10], prices[10], marker='o', color='red', markersize=10)
        ax1.set_title('Price vs K')
        ax1.set_xlabel('K')
        ax1.set_ylabel('Price')
        ax1.grid(True)
        st.pyplot(fig1)





        # Trinomial tree visualization
        st.subheader('Trinomial Tree Visualization')
        fig2 = plot_trinomial_tree(r)
        st.pyplot(fig2)