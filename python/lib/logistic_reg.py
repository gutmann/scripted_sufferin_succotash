import statsmodels.api as sm
import numpy as np 

class LogitWeight(sm.Logit):
    def __init__(self, endog, exog, **kargs):
        super(sm.Logit, self).__init__(endog, exog, **kargs)
        self.weights = kargs.get("weights", np.ones_like(endog))

    def loglike(self, params):
        q = 2*self.endog - 1
        X = self.exog
        return np.sum(self.weights*np.log(self.cdf(q*np.dot(X, params))))

    def loglikeobs(self, params):
        q = 2*self.endog - 1
        X = self.exog
        return self.weights*np.log(self.cdf(q*np.dot(X, params)))

    def jac(self, params):
        y = self.endog
        X = self.exog
        L = self.cdf(np.dot(X, params))
        return ((y - L) * self.weights)[:, None] * X

    def score(self, params):
        y = self.endog
        X = self.exog
        L = self.cdf(np.dot(X, params))
        return np.dot((y - L)*self.weights, X)

    def hessian(self, params):
        X = self.exog
        L = self.cdf(np.dot(X, params))
        return -np.dot((self.weights*L*(1-L)*X.T), X)
