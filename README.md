# Heston-Model-
def function (S0, mu, v0, rho, kappa, theta, etha, T, dt) :

	# Generate a Montcarlo simulation for the Heston Model

	# Generate random Brownian Motion

	Mu = np.array([0, 0])
	COV = np.matrix([[1, rho], [rho, 1]])
	W = np.random.multivariate_normal(MU, COV, T)
        W_S = W[:,0]
        W_v = W[:,1]

        # Generate paths
    vt    = np.zeros(T)
    vt[0] = v0
    St    = np.zeros(T)
    St[0] = S0
    for t in range(1,T):
        vt[t] = np.abs(vt[t-1] + kappa*(theta-(vt[t-1]))*dt + etha*np.sqrt((vt[t-1]))*W_v[t]*np.sqrt(dt))
        St[t] = St[t-1]*np.exp((mu - 0.5*vt[t-1])*dt + np.sqrt(vt[t-1]*dt)*W_S[t])

    return St, vt

