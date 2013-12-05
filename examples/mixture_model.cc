#include "chainmcmc/chainmcmc.hh"

using namespace chainmcmc;

/* Here we add the mixture model discussed in section 3.2 of:
 *
 		Carlin, Bradley P., and Siddhartha Chib. 1995. “Bayesian Model Choice via Markov Chain Monte Carlo Methods.” Journal of the Royal Statistical Society. Series B (Methodological) 57 (3) (January 1): 473–484.
 * 
 * This allows us to test our results versus previous results and make sure we 
 * didn't make any mistakes.
 * Table 3. in that artikel shows their findings
 */

int main() {
	// Data
	std::vector<double> velocities = { 
		9172, 9350, 9483, 9558, 9775, 10227,
		10406, 16084, 16170, 18419, 18552, 18600,
		18927, 19052, 19070, 19330, 19343, 19349,
		19440, 19473, 19529, 19541, 19547, 19663,
		19846, 19856, 19863, 19914, 19918, 19973,
		19989, 20166, 20175, 20179, 20196, 20215,
		20221, 20415, 20629, 20795, 20821, 20846,
		20875, 20986, 21137, 21492, 21701, 21814,
		21921, 21960, 22185, 22209, 22242, 22249,
		22314, 22374, 22495, 22746, 22747, 22888,
		22914, 23206, 23241, 23263, 23484, 23538,
		23542, 23666, 23706, 23711, 24129, 24285,
		24289, 24366, 24717, 24990, 25633, 26960,
		26995, 32065, 32789, 34279 };
	
	// Setup likelihood (can probably just depend on length of pars to work out number of mixtures)
	likelihood_t ll = [&velocities]( const std::vector<parameter_t> &pars )  {
		// pars consists of parameters for each mixture, 
		// with each mixture being represented by mu, weight
		// First parameter gives the sigma for all mixtures
		//
		// The last weight is 1-sum other weights
		std::vector<prior_t> mixtures;
		double sum_weight = 0;
		for ( size_t i = 1; i < velocities.size(); i += 2 ) {
			mixtures.push_back( [&pars, &i, &sum_weight]( const parameter_t &v ) {
					if (i+1==pars.size())
						return (1.0-sum_weight)*prior::normal( pars[i], pars[0] )( v );
					else {
						sum_weight += pars[i+1];
						return pars[i+1]*prior::normal( pars[i], pars[0] )( v );
					}
			} );
		}
		double myll = 0;
		for ( auto &v : velocities ) {
			double sum = 0;
			for ( auto &n : mixtures ) {
				sum += n( v );
			}
			myll += log( sum );
		}

		return myll;
	};

	return 0;
}
