from __future__ import division
import numpy as np


class Routing(object):
    """
    A hydrologic engine for a water management and decision support system.

    This is the central class of the hydrologic component. It is initialized with the
    connectivity matrix of the network and routes runoff generated at each node subcatchment
    through out the network.

    It uses outputs from the rainfall runoff models as water so

    """
    def __init__(self, conn, dt):
        """
        Class initialization. It requires a connectivity matrix describing the topology of the
         water distribution network. The connectivity matrix is square and sparse, of size
           equal to the number of nodes to the network. Each element in the matrix represents
           a link (reach, channel) between nodes (arbitrary points in the river network,
           diversion points, gauges, water users, etc).

        :param conn: nxn binary connectivity matrix, where n is the number of nodes in the
        network.
        :param dt: time step (seconds, scalar).

        """

        self.conn = conn
        self.dt = dt
        self.n = np.shape(conn)[0]

    def muskingum_routing(self, Qt, K, e, qnew, qold):
        """Routes water through a network graph using the Muskingum-Cunge method.

        This function takes an initial distribution of streamflows at nodes
        in the network and routes them one time step.

        :param Qt: vector of streamflows at time t
        :param K: vector of reach storage parameters of size n.
        :param e: vector of balances between inflows and outflows
        :param q: vector of lateral inflows

        :return: vector of size n with streamflows at time t+1.
        """

        # Add here code to check the stability of the M-C algorithm
        # Only the upper bound this condition is necessary for stability,
        # to lower bound for dt is to avoid possible negative flows
        dt = self.dt
        max_stab = 2 * K * (1 - e)
        min_stab = 2 * K * e  # implementation of this condition may require that K is halved

        # Adjust dt to maintain stability of Muskingum solution
        # The algorithm halves the time step until it meets the criterium for all reaches
        # so the shortest/fastest reach will control the integration time step
        n = 1 # maintains the number of fractions that need ot be computed to simulate the full day
        while np.any(np.greater(self.dt, max_stab)):
            dt /= 2
            n *= 2

        # adjust Q and q to the reduced dt
        Qt = Qt/n
        qnew = qnew/n
        qold = qold/n

        # while np.any(np.greater(min_stab, self.dt)):
        #     K = K/2 # assuming K = dx/c
        #     # e = do something about e here, which is also a function of dx
        #     n *= 2
        #     min_stab = 2 * K * e

        for i in range(n):

            a = np.diag(K * (1-e) + dt * 0.5)
            b = np.diag(K * e - dt * 0.5)
            c = np.diag(K * (1-e) - dt * 0.5)
            d = np.diag(K * e + dt * 0.5)

            lhs = (a+np.dot(self.conn, b)).T
            rhs = np.dot((d+np.dot(self.conn, c)).T, Qt) + np.diag(a)*qnew#+ np.diag(d)*qold - np.diag(b)*qnew

           # print lhs
           # print rhs

            Qt1 = np.linalg.solve(lhs, rhs)
            Qt = Qt1

            #print(np.allclose(np.dot(lhs, Qt1), rhs))

        return Qt1




