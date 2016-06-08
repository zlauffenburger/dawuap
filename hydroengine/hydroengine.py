from __future__ import division
import numpy as np

class routing(object):
    """A rainfall-runoff model for water management and decision support system

    """
    def __init__(self, conn, dt):
        """Cass hydro.

        :args:
            :conn: nxn connectivity matrix
            :dt: time step (scalar).

        """

        self.conn = conn
        self.dt = dt
        self.n = np.shape(conn)[0]



    def muskingum_routing(self, Qt, K, e, q):
        """Routes water through a network graph using the Muskingum-Cunge method.

        This function takes an initial distribution of streamflows at noddes
        in the network and routes them one time step.

        :args:
            :Qt: vector of streamflows at time t
            :K: vector of reach storage parameters of size n.
            :e: vector of balances between inflows and outflows
            :q: vector of lateral inflows

        :return: vector of size n with streamflows at time t+1.
        """

        a = np.diag(K*(1-e) + self.dt*0.5)
        b = np.diag(K*e-self.dt*0.5)
        c = np.diag(K*(1-e) - self.dt*0.5)
        d = np.diag(K * e - self.dt * 0.5)

        lhs = (a+np.dot(self.conn,b)).T
        rhs = np.dot((d+np.dot(self.conn,c)).T,Qt) + qstuff

        print lhs
        print rhs

        Qt1 = np.linalg.solve(lhs,rhs)

        print(np.allclose(np.dot(lhs,Qt1), rhs))

        return Qt1




