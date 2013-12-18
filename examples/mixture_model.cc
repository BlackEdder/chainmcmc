#include "chainmcmc/chainmcmc.hh"

#include<realtimeplot/realtimeplot.h>

using namespace chainmcmc;

/* Here we add the mixture model discussed in section 3.2 of:
 *
 		Carlin, Bradley P., and Siddhartha Chib. 1995. “Bayesian Model Choice via Markov Chain Monte Carlo Methods.” Journal of the Royal Statistical Society. Series B (Methodological) 57 (3) (January 1): 473–484.
 * 
 * This allows us to test our results versus previous results and make sure we 
 * didn't make any mistakes.
 * Table 3. in that artikel shows their findings
 */

joint_prior_t convert_to_joint_prior( const std::vector<double> & init ) {
	// pars consists of parameters for each mixture, 
	// with each mixture being represented by mu, weight
	// First parameter gives the sigma for all mixtures
	//
	// The last weight is 1-sum other weights
	

	double alpha = 2+pow( init[1], 2 )/pow( init[2], 2 );
	double beta = init[1]+pow( init[1], 3 )/pow( init[2], 2 );
	prior_t pr_var = prior::inverse_gamma( alpha, beta );

	std::vector<prior_t> pr_means;
	for ( size_t i = 0; i < init[0]; ++i ) {
		pr_means.push_back( prior::normal( init[ (i*2)+3 ], init[ (i*2)+4 ] ) );
	}

	// Dir for the weights
	std::vector<double> alphas( init[0] );
	double sum_alphas = 1;
	alphas[0] = init[ 3+2*init[0] ]*sum_alphas;
	size_t j = 1;
	for ( size_t i = 5+2*init[0]; i < init.size(); i += 2 ) {
		alphas[j] = init[i]*sum_alphas;
		++j;
	}

	// scale alphas to get the right variance
	// alphas[0] *= sa.. since alphas[0] == 1, alphas[0] = sa
	// var x0 = ab0*sa*(sa*sum_alphas-ab0*sa)/(pow(sa*sum_alphas,2)*(sa*sum_alphas+1));
	// sd x0 = init[4+2*init[0]]
	double sa = - (pow(sum_alphas,2)*pow(init[4+2*init[0]],2)-alphas[0]*sum_alphas+pow(alphas[0],2))/(pow(sum_alphas,3)*pow(init[4+2*init[0]],2));
	for ( auto &alp : alphas )
		alp *= sa;

	joint_prior_t pr_dir = prior::dirichlet( alphas );
	size_t dim = init[0];

	joint_prior_t jp( [pr_var, pr_means, pr_dir, dim ] ( 
				const std::vector<parameter_t> &pars ) {
		double pr = pr_var( pow(pars[0],2) );

		std::vector<parameter_t> par_alphas;
		for ( size_t i = 0; i < dim; ++i ) {
			pr *= pr_means[i]( pars[(i*2)+1] ); 

			if ((i*2)+2<pars.size()) {
				if (pars[(i*2)+1]>pars[(i*2)+3]) // Keep mixture means ordered/sorted
					return 0.0;
				par_alphas.push_back( pars[ (i*2)+2 ] );
			}
		}

		pr *= pr_dir( par_alphas );


		return pr;
	} );
	return jp;
}

double dmixture( const double &v, const std::vector<parameter_t> &pars ) {
	std::vector<prior_t> mixtures;
	double sum_weight = 0;
	for ( size_t i = 1; i < pars.size(); i += 2 ) {
		if (i+1<pars.size())
			sum_weight += pars[i+1];
		mixtures.push_back( [&pars, i, sum_weight]( const parameter_t &v ) {
				if (i+1<pars.size())
				return pars[i+1]*prior::normal( pars[i], pars[0] )( v );
				else 
				return (1.0-sum_weight)*prior::normal( pars[i], pars[0] )( v );
				} );
	}
	double sum = 0;
	for ( auto &n : mixtures ) {
		sum += n( v );
	}
	return sum;
}

void plot_results( const std::vector<parameter_t> &pars ) {
	realtimeplot::PlotConfig conf;
	conf.fixed_plot_area = true;
	conf.min_x = 0;
	conf.max_x = 40000;
	conf.max_y = 2e-4;
	conf.aspect_ratio = 2.5;
	realtimeplot::Plot pl = realtimeplot::Plot(conf);
	for ( size_t i = 0; i < 40000; i += 100 ) {
		pl.line_add( i, dmixture( i, pars ) );
	}
}

double accept_prob( const std::vector<parameter_t> &new_pars,
		const std::vector<parameter_t> &old_pars, const likelihood_t &ll, 
		const joint_prior_t &jp ) {
	double alpha = exp( step::mh_log_weight( ll, new_pars, jp, 1 ) -
			step::mh_log_weight( ll, old_pars, jp, 1 ) );
	//std::cout << alpha << std::endl;
	alpha = std::min( 1.0, alpha );
	//std::cout << "BLA: " << alpha << std::endl;
	return alpha;
}

double log_marginal_likelihood_vit( const likelihood_t &ll, 
		const joint_prior_t &jp, 
		const std::vector<std::vector<parameter_t> > &tr ) {
	std::mt19937 eng;
	auto theta_star = trace::means( tr );
	auto theta_var = trace::variances_sample( tr );
	size_t no_samples = 2000;
	double numer = 0;
	double denom = 0;
	std::cout << "Star: " << theta_star << std::endl;
	std::cout << "Var:  " << theta_var << std::endl;
	for (size_t i = 0; i < theta_star.size(); ++i) {
		theta_star[i] -= sqrt(theta_var[i]);
	}
	std::cout << "Star: " << theta_star << std::endl;

	for (size_t i = 0; i < no_samples; ++i) {
		double prod_numer = 1;
		double prod_denom = 1;
		for (size_t j = 0; j < theta_star.size(); ++j) {
			auto theta_g = theta_star;
			theta_g[j] = tr[ tr.size()-1-i ][j];
			prod_numer *= prior::normal(theta_star[j], (theta_var[j]))(theta_g[j]);
			prod_numer *= accept_prob( theta_star, theta_g, ll, jp );
			double sum_denom = 0;
			for (size_t m = 0; m < 100; ++m) {
				auto theta_g = theta_star;
				std::normal_distribution<double> rnorm( 0, sqrt(theta_var[j]) );
				theta_g[j] += rnorm( eng );
				sum_denom += accept_prob( theta_g, theta_star, ll, jp );
			}
			prod_denom *= sum_denom/1000.0;
		}
		numer += prod_numer;
		denom += prod_denom;
	}
	std::cout << numer << " / " << denom << std::endl;
	return step::mh_log_weight( ll, theta_star, jp, 1 )-log(1.0/no_samples*numer/denom);
}

double log_marginal_likelihood( std::mt19937 &eng, const likelihood_t &ll, 
		const std::vector<parameter_t> &init_pars, const 
		joint_prior_t &jp ) {
	std::vector<double> as = { 0,0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 
		0.8, 0.9,	1.0 };

	double exp_ll = 0;
	double sum_ll = 0;

	for ( size_t i = 0; i < as.size(); ++i ) {
		double ti = pow(as[i],5);
		likelihood_t hot = [&ll, &ti]( const std::vector<parameter_t> &pars ) {
			if (ti == 0)
				return 0.0;
			else {
				return ti*ll(pars);
			}
		};
		auto chain = spawn<Chain>( eng, hot, init_pars, jp );

		std::stringstream out;
		auto logger = spawn<Logger>( out );
		send( chain, atom("logger"), logger );

		send( chain, atom("run"), 100000, false );
		send( chain, atom("no_adapt") );
		send( chain, atom("run"), 100000, true );
		send( chain, atom("log_weight") );
		receive( 
				on( arg_match ) >> [&out, &ti]( const double &lw ) {
				//std::cout << out.str() << std::endl;
				std::cout << ti << std::endl;
				std::cout << "Weight: " << exp( lw ) << std::endl; } );

		auto tr = trace::read_trace_per_sample( out );

		double mean_ll = 0;
		size_t dim = tr.size();
		for ( auto & sample : tr ) {
			auto l = ll( sample );
			if (std::isfinite(l))
				mean_ll += l;
			else
				--dim;
		}
		mean_ll /= dim;
		if (!std::isfinite( mean_ll))
			mean_ll = (-std::numeric_limits<double>::max())/100.0;
		std::cout << "Mean ll " << mean_ll << std::endl;

		if ( i != 0 ) {
			sum_ll += (pow(as[i],5)-pow(as[i-1],5))*(mean_ll+exp_ll)/2.0;
		}
		exp_ll = mean_ll;
	}

	return sum_ll;
}

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
	
	// Prior setup, taken from the article. First par is number of mixtures.
	// Then the mean and sd for the prior of the variance
	// Then for each mixture the priors for different means defined by a mean
	// and a sd.
	// Finally the means and sds for the weights (number of weights - 1)
	std::vector<double> init_prior1 = { 3, 20000, 20000, 
		9000, 5000, 18000, 5000, 30000, 5000, 
		1.0/3, 0.236, 1.0/3, 0.236, 1.0/3, 0.236 };
	std::vector<double> init_prior2 = { 4, 15000, 15000, 
		9000, 5000, 18000, 5000, 22000, 5000, 30000, 5000, 
		0.136, 0.100, 0.364, 0.139, 0.364, 0.139, 0.136, 0.100 };

	joint_prior_t jp1 = convert_to_joint_prior( init_prior1 );
	joint_prior_t jp2 = convert_to_joint_prior( init_prior2 );

	std::vector<parameter_t> init_pars1 = { 20000, 9000, 1.0/3, 18000, 1.0/3, 30000 };
	std::vector<parameter_t> init_pars2 = { 15000, 9000, 1.0/4, 18000, 1.0/3, 22000, 1.0/3, 30000 };
	
	// Setup likelihood (can probably just depend on length of pars to work out number of mixtures)
	likelihood_t ll = [&velocities]( const std::vector<parameter_t> &pars )  {
		// pars consists of parameters for each mixture, 
		// with each mixture being represented by mu, weight
		// First parameter gives the sigma for all mixtures
		//
		// The last weight is 1-sum other weights
		std::vector<prior_t> mixtures;
		double sum_weight = 0;
		for ( size_t i = 1; i < pars.size(); i += 2 ) {
			if (i+1<pars.size())
				sum_weight += pars[i+1];
			mixtures.push_back( [&pars, i, sum_weight]( const parameter_t &v ) {
					if (i+1<pars.size())
						return pars[i+1]*prior::normal( pars[i], pars[0] )( v );
					else 
						return (1.0-sum_weight)*prior::normal( pars[i], pars[0] )( v );
			} );
		}
		double myll = 0;
		for ( auto &v : velocities ) {
			double sum = 0;
			for ( auto &n : mixtures ) {
				sum += n( v );
			}
			//std::cout << v << " " << sum << " " << myll << std::endl;
			if ( sum == 0 )
				myll -= std::numeric_limits<double>::max()/1000.0;
			else 
				myll += log( sum );
		}
		return myll;
	};

	std::mt19937 eng;

	std::cout << "Friel:" << std::endl;
	auto ml1 = log_marginal_likelihood( eng, ll, init_pars1, jp1 );
	auto ml2 = log_marginal_likelihood( eng, ll, init_pars2, jp2 );
	std::cout << ml1 << std::endl;
	std::cout << ml2 << std::endl;;
	std::cout << ml2-ml1 << " " << ml1-ml2 << std::endl;
	std::cout << exp(ml2-ml1) << " " << exp(ml1-ml2) << std::endl;
	std::cout << "End Friel:" << std::endl;



	std::stringstream out;
	auto chainC = ChainController( ll, init_pars1, jp1, 100000, 300000, 8, out );

	auto tr1 = trace::read_trace_per_sample( out );
	auto means1 = trace::means( tr1 );
	std::cout << std::endl << std::endl << step::mh_log_weight( ll, means1, jp1, 1 ) << 
		std::endl << ll(means1) << ": " << means1 << std::endl;
	plot_results( means1 );
	auto s_vars = trace::variances_sample( tr1 );
	std::transform( s_vars.begin(), s_vars.end(), s_vars.begin(),
			[]( const double &el ) {
			return sqrt(el); } );
	std::cout << s_vars << std::endl;
	std::cout << std::endl << std::endl;

	out.clear();
	chainC = ChainController( ll, init_pars2, jp2, 100000, 300000, 8, out );

	auto tr2 = trace::read_trace_per_sample( out );
	auto means2 = trace::means( tr2 );
	std::cout << std::endl << std::endl << step::mh_log_weight( ll, means2, jp2, 1 ) << 
		std::endl << ll(means2) << ": " << means2 << std::endl;
	plot_results( means2 );
	s_vars = trace::variances_sample( tr2 );
	std::transform( s_vars.begin(), s_vars.end(), s_vars.begin(),
			[]( const double &el ) {
			return sqrt(el); } );
	std::cout << s_vars << std::endl;


	ml1 = log_marginal_likelihood_vit( ll, jp1, tr1 );
	ml2 = log_marginal_likelihood_vit( ll, jp2, tr2 );

	std::cout << ml1 << std::endl;
	std::cout << ml2 << std::endl;;
	std::cout << ml2-ml1 << " " << ml1-ml2 << std::endl;
	std::cout << exp((ml1+1-ml2)/2) << " " << exp(-ml2) << std::endl;





	return 0;
}
