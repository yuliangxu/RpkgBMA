
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]


#include "BMA_header.h"

// -----------------------GS update to give initial values--------------------------------- //
//' Initial values for image-on-scalar regression
//' 
//' @param data an R list object that contains the data
//' @param init R list object with initial values
//' @param region_idx_cpp a list of region idx for all regions starting from 0
//' @param kernel a list object containing Q and D
//' @param n_mcmc number of iterations
//' @param display_progress bool, 1 = display progress bar
// [[Rcpp::export]]
Rcpp::List M_regression_GS(Rcpp::List& data, Rcpp::List& init,
                     Rcpp::List& region_idx_cpp,
                     Rcpp::List& kernel, arma::uword n_mcmc,
                     bool display_progress = true){
  // input data
  
  Rcpp::Timer timer;  
  timer.step("start of precomputation"); 
  arma::colvec X = data["X"];
  arma::mat C = data["C"]; 
  arma::mat M = data["M"];
  arma::uword n = X.n_elem;
  // input initial parameters
  arma::colvec zetam = init["zetam"];
  double sigma_alpha = init["sigma_alpha"]; 
  double sigma_M = init["sigma_M"];
  arma::mat theta_eta = init["theta_eta"]; 
  double sigma_eta = init["sigma_eta"];
  arma::colvec theta_alpha = init["theta_alpha"];
  double a = init["a"]; double b = init["b"];
  double sigma_zetam = init["sigma_zetam"];
  double sigma_alpha2_inv = 1/sigma_alpha/sigma_alpha;
  double sigma_M2_inv = 1/sigma_M/sigma_M;
  double sigma_zetam2_inv = 1/sigma_zetam/sigma_zetam;
  double sigma_eta2_inv = 1/sigma_eta/sigma_eta;
  double X2_sum = dot(X,X);
  
  // input basis functions
  Rcpp::List D = kernel["D"];Rcpp::List Q = kernel["Q"];
  
  
  // get M_star, D_vec, p, L from basis
  arma::uword num_region = region_idx_cpp.length();
  arma::uvec p_length; arma::uvec L_all;
  p_length.set_size(num_region); L_all.set_size(num_region);
  arma::mat M_star;
  arma::colvec D_vec; arma::colvec q_vec;
  for(arma::uword i=0; i<num_region; i++){
    arma::mat Q_i = Q[i];
    p_length(i) = Q_i.n_rows; L_all(i)=Q_i.n_cols;
    arma::uvec  idx = region_idx_cpp[i]; 
    M_star = arma::join_cols(M_star, Q_i.t()*M.rows(idx));
    arma::colvec D_i = D[i];
    D_vec = arma::join_cols(D_vec, D_i);
    arma::rowvec q_i = sum(Q_i,0);
    q_vec = arma::join_cols(q_vec, q_i.t());
  }
  arma::uword L = sum(L_all); arma::uword p = sum(p_length);
  double q2_sum = dot(q_vec,q_vec);
  // M_star = M_star/sqrt(p);
  arma::mat M_star_t = M_star.t();
  
  arma::uword m = zetam.n_elem;
  arma::mat C_t = C.t();
  
  // initialize mcmc sequences
  arma::mat theta_alpha_mcmc = arma::zeros(L,n_mcmc );
  arma::mat zetam_mcmc = arma::zeros(m,n_mcmc );
  arma::colvec sigma_M2_inv_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_alpha2_inv_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_zetam2_inv_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_eta2_inv_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec logLL_mcmc = arma::zeros(n_mcmc,1);
  arma::mat zeta_term =  q_vec*zetam.t()*C;
  Progress prog(n_mcmc*num_region, display_progress);
  timer.step("start of iteration"); 
  for(arma::uword iter=0; iter< n_mcmc; iter++){
    prog.increment(); 
    
    // 1. update theta_alpha
    arma::mat M_res = M_star - zeta_term - theta_eta;
    arma::colvec Sigma_theta_alpha = 1/(sigma_alpha2_inv/D_vec + sigma_M2_inv*X2_sum);
    M_res.each_row() %= X.t();
    arma::colvec mu_theta_alpha = sum(M_res,1)%Sigma_theta_alpha*sigma_M2_inv;
    theta_alpha =  arma::randn(L,1)%sqrt(Sigma_theta_alpha) + mu_theta_alpha;
    theta_alpha_mcmc.col(iter) =  mu_theta_alpha;

    // 2. zetam

    arma::mat Sigma_zetam_inv = sigma_zetam2_inv + sigma_M2_inv*q2_sum*C*C_t;
    arma::mat Sigma_zetam = inv_sympd(Sigma_zetam_inv );
    arma::mat zeta_res = M_star - theta_alpha*X.t() - theta_eta;
    arma::rowvec temp_zeta = q_vec.t()*zeta_res;
    arma::rowvec mu_zetam = temp_zeta*C_t *Sigma_zetam* sigma_M2_inv;
    zetam = arma::mvnrnd( mu_zetam.t(), Sigma_zetam);
    zetam_mcmc.col(iter) = zetam;


    // 3. theta_eta
    zeta_term =  q_vec*zetam.t()*C;
    arma::mat eta_res = M_star - theta_alpha*X.t() - zeta_term; //L x n
    arma::colvec eta_Sigma_vec = 1/(sigma_M2_inv + sigma_eta2_inv/D_vec);
    arma::colvec mu_eta_i;
    // mat mu_eta;
    for(arma::uword i=0;i<n;i++){
      mu_eta_i = (eta_res.col(i)*sigma_eta2_inv)%eta_Sigma_vec;
      arma::colvec theta_eta_i = arma::randn(L)%arma::sqrt(eta_Sigma_vec)+ mu_eta_i;
      theta_eta.col(i) = theta_eta_i - arma::mean(theta_eta_i);
      // mu_eta = join_rows(mu_eta,mu_eta_i);
      // theta_eta.col(i) = mu_eta;
    }


    // 4. sigma_alpha
    sigma_alpha2_inv = arma::randg( arma::distr_param(a + L*0.5, 1/(b + dot(theta_alpha,theta_alpha/D_vec)/2)) );
    sigma_alpha2_inv_mcmc(iter) = sigma_alpha2_inv;
    // 5. sigma_M
    arma::mat M_reg_res = M_star - theta_alpha*X.t() - zeta_term - theta_eta;
    double M_norm = norm(M_reg_res,"fro");
    sigma_M2_inv = arma::randg( arma::distr_param(a + n*L/2,1/(b + M_norm*M_norm/2) ) );
    sigma_M2_inv_mcmc(iter) = sigma_M2_inv;

    //   6. sigma_zetam -> half-cauchy
    // sigma_zetam2_inv = arma::randg( arma::distr_param(a + zetam.n_elem/2,1/(b + dot(zetam,zetam)/2) ) );
    // sigma_zetam2_inv_mcmc(iter) = sigma_zetam2_inv;

    // //   7. sigma_eta -> giving too large variance
    double eta_norm = norm(theta_eta,"fro");
    sigma_eta2_inv = arma::randg( arma::distr_param(a + 0.5*n*L,1/(b + eta_norm*eta_norm/2)) );
    sigma_eta2_inv_mcmc(iter) = sigma_eta2_inv;
    
    // logLL
    arma::mat res = M_star - zeta_term - theta_eta - theta_alpha*X.t();
    double res_norm = norm(M_star,"fro");
    logLL_mcmc(iter) = -0.5*sigma_M2_inv*res_norm*res_norm;
  }
  timer.step("end of iterations");     
  
  return Rcpp::List::create(Rcpp::Named("theta_alpha_mcmc")=theta_alpha_mcmc,
                            Rcpp::Named("theta_eta") = theta_eta,
                              Rcpp::Named("zetam_mcmc")= zetam_mcmc,
                              Rcpp::Named("sigma_M2_inv_mcmc")= sigma_M2_inv_mcmc ,
                              Rcpp::Named("sigma_alpha2_inv_mcmc")= sigma_alpha2_inv_mcmc,
                              Rcpp::Named("sigma_zetam2_inv_mcmc")=  sigma_zetam2_inv_mcmc,
                              Rcpp::Named("sigma_eta2_inv_mcmc")= sigma_eta2_inv_mcmc,
                              Rcpp::Named("logLL_mcmc")= logLL_mcmc,
                              Rcpp::Named("Timer")=timer);
}

// // -----------------------block update theta_alpha--------------------------------- //
// [[Rcpp::export]]
Rcpp::List M_regression_region_block(arma::colvec& Y, arma::mat& M,
                               arma::colvec& X, arma::mat& C, arma::uvec L_all,
                               arma::uword num_region, Rcpp::List& region_idx,
                               int n_mcmc, Rcpp::List& K, int stop_burnin,
                               int unadjusted_langevin,
                               double lambda, double& target_accept,
                               Rcpp::List& init,
                               int interval,
                               double step,
                               bool display_progress = true){
  Rcpp::Timer timer;
  BMA bma;
  bma.set_seed(1);
  timer.step("start of precomputation");
  arma::uword p = M.n_rows;
  arma::uword n = M.n_cols;
  
  arma::uvec L_cumsum = cumsum(L_all);
  arma::uword L_max = L_cumsum(num_region-1);
  // input
  arma::colvec theta_alpha = init["theta_alpha"];
  arma::colvec D = init["D"];
  arma::colvec zetam = init["zetam"];
  arma::uword m = zetam.n_elem;
  double sigma_M = init["sigma_M"], sigma_M2 = sigma_M*sigma_M, sigma_M2_inv = 1/sigma_M2;
  double sigma_alpha = init["sigma_alpha"], sigma_alpha2 = sigma_alpha*sigma_alpha, sigma_alpha2_inv = 1/sigma_alpha2;
  double sigma_eta = init["sigma_eta"], sigma_eta2_inv = 1/sigma_eta/sigma_eta;
  arma::mat theta_eta = init["theta_eta"];double sigma_zetam2_inv = 100;
  arma::colvec step_all = step*arma::ones(num_region);
  double sigma_lambda = 0.5;
  // hyper parameter for inverse-Gamma
  double a = 1.0, b=1.0;
  if(C.n_cols !=zetam.n_rows){
    Rcpp::Rcout<<"Error: dimensions of C and zetam don't match!"<<
      "dim of C = "<<size(C)<<"; dim of zetam = "<<size(zetam)<<std::endl;
    return Rcpp::List::create(Rcpp::Named("ERROR")=1);
  }
  arma::mat C_t = C.t();
  Rcpp::List M_star_pre_eta(num_region);
  arma::colvec alpha = arma::zeros(p,1);
  arma::colvec q_vec;
  for(int l=0; l<num_region; l++){
    arma::uvec idx = region_idx[l];
    arma::mat Q = K[l]; 
    // M_star_pre_eta[l] = Q_t*(M.rows(idx) - arma::ones(idx.n_elem,1)*zetam.t()*C);
    arma::uvec L_idx;
    if(l==0){
      L_idx = arma::linspace< arma::uvec>(0,L_cumsum(l)-1,L_all(l));
    }else{
      L_idx = arma::linspace< arma::uvec>(L_cumsum(l-1),L_cumsum(l)-1,L_all(l));
    }
    alpha(idx) = Q*theta_alpha(L_idx);
    arma::rowvec q_i = sum(Q,0);
    q_vec = arma::join_cols(q_vec, q_i.t());
  }
  double q2_sum = dot(q_vec,q_vec);
  //return
  arma::colvec lambda_mcmc = arma::zeros(n_mcmc);
  arma::mat theta_alpha_mcmc = arma::zeros(L_max,n_mcmc);
  arma::colvec logll_mcmc = arma::zeros(n_mcmc);
  arma::mat track_grad_f = arma::zeros(L_max,n_mcmc);
  arma::mat track_rho = arma::zeros(num_region, n_mcmc);
  arma::mat track_step = arma::zeros(num_region,n_mcmc);
  arma::colvec accept_lambda = arma::zeros(n_mcmc);
  arma::cube track_rho_compo = arma::zeros(num_region,4, n_mcmc);
  arma::mat emp_accept = arma::zeros( n_mcmc/interval,num_region);
  arma::colvec accept = arma::zeros(n_mcmc*num_region);
  arma::mat accept_block = arma::zeros(n_mcmc,num_region);
  arma::uword all_iter=0;
  arma::uword num_block = num_region;
  
  // GS: initialize mcmc sequences
  arma::mat zetam_mcmc = arma::zeros(m,n_mcmc );
  arma::colvec sigma_M2_inv_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_alpha2_inv_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_zetam2_inv_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_eta2_inv_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec logLL_mcmc = arma::zeros(n_mcmc,1);
  arma::mat zeta_term =  q_vec*zetam.t()*C_t;
  
  Progress prog(n_mcmc*num_region, display_progress);
  timer.step("start of iteration");
  for(int iter=0; iter<n_mcmc; iter++){
      double logll_M = 0;
  
      if(iter==stop_burnin){
        // start the timer
        timer.step("stop of burnin");        // record the starting point
      }
  //     // Rcout<<"test 3"<<std::endl;
  //     
  arma::colvec grad_f_all = arma::zeros(L_max,1);
      arma::colvec rho_all = arma::zeros(num_block,1);
      arma::mat rho_compo = arma::zeros(num_block,4);
  
      // check if region_idx starts from 0!
      // start block update
      for(arma::uword m=0; m < num_region; m++){
        prog.increment();
        // Rcout<<"iter="<<iter<<"; region="<<m<<std::endl;
        arma::uvec delta = find(abs(alpha)>lambda);
        arma::uvec idx = region_idx[m];
        arma::uvec delta_Q = find(abs(alpha(idx))>lambda);
        // Rcout<<"test 4"<<std::endl;
        arma::mat Q = K[m];
        arma::uvec L_idx;
        if(m==0){
          L_idx = arma::linspace< arma::uvec>(0,L_cumsum(m)-1,L_all(m));
        }else{
          L_idx = arma::linspace< arma::uvec>(L_cumsum(m-1),L_cumsum(m)-1,L_all(m));
        }
  //       // Rcout<<"test 5"<<std::endl;
        arma::uvec delta_in_block = intersect(idx,delta);

        arma::mat K_block_t = Q.t();
        arma::mat M_star = K_block_t*(M.rows(idx) - arma::ones(idx.n_elem,1)*zetam.t()*C_t) - theta_eta.rows(L_idx);
        arma::mat temp = M_star - K_block_t.cols(delta_Q)*(alpha.rows(delta_in_block)-sign(alpha.rows(delta_in_block))*lambda)*X.t();

        arma::mat temp_X = temp;
        temp_X.each_row()%= X.t();
        arma::colvec temp_sum = arma::sum(temp_X,1)*sigma_M2_inv;
        arma::colvec grad_f = -theta_alpha(L_idx)/D(L_idx)*sigma_alpha2_inv+
          K_block_t.cols(delta_Q)*Q.rows(delta_Q) *temp_sum; // L x 1, use smooth function in grad

        double step = step_all(m);
        arma::colvec theta_alpha_diff = step*grad_f+sqrt(2*step)*arma::randn(size(grad_f));
        arma::colvec theta_alpha_new_block = theta_alpha(L_idx) + theta_alpha_diff;
        grad_f_all(L_idx) = grad_f;

        if(iter< unadjusted_langevin){
          accept(all_iter) = 1;
          accept_block(iter,m) = 1;
          alpha(idx) = alpha(idx) + Q*theta_alpha_diff;
          theta_alpha(L_idx) = theta_alpha_new_block;

        }else{
          // MH step
          // Rcout<<"test 6"<<std::endl;
          double log_target_density = -0.5*bma.square(norm(theta_alpha(L_idx)/sqrt(D(L_idx)),2))*sigma_alpha2_inv-
           0.5*bma.square(arma::norm(temp,"fro"))*sigma_M2_inv;


          arma::colvec alpha_new = alpha;
          alpha_new(idx) = alpha_new(idx) + Q*theta_alpha_diff;
          arma::uvec delta_new = find(abs(alpha_new)>lambda);
          arma::uvec delta_Q_new = find(abs(alpha_new(idx))>lambda);
          arma::uvec delta_in_block_new = intersect(idx,delta_new);
          arma::mat temp_new = M_star - K_block_t.cols(delta_Q_new)*(alpha_new.rows(delta_in_block_new)-sign(alpha_new.rows(delta_in_block_new))*lambda)*X.t();

          arma::mat temp_X_new = temp_new;
          temp_X_new.each_row()%= X.t();
          arma::colvec temp_sum_new = arma::sum(temp_X_new,1)*sigma_M2_inv;
          arma::colvec grad_f_new = -theta_alpha_new_block/D(L_idx)*sigma_alpha2_inv+
            K_block_t.cols(delta_Q_new)*Q.rows(delta_Q_new) *temp_sum_new; // L x 1, use smooth function in grad
          // Rcout<<"test 6-5"<<std::endl;
          double log_target_density_new = -0.5*bma.square(norm( theta_alpha_new_block/sqrt(D(L_idx)),2))*sigma_alpha2_inv-
            0.5*bma.square(arma::norm(temp_new,"fro"))*sigma_M2_inv;
  //         // Rcout<<"test 7"<<std::endl;
          double log_q = -1/4/step * bma.square(norm(-theta_alpha_diff-step*grad_f_new,2));
          double log_q_new = -1/4/step * bma.square(norm(theta_alpha_diff-step*grad_f,2));
          double rho = log_target_density_new + log_q - log_target_density - log_q_new;
          rho_all(m) = rho;
          rho_compo.row(m) = {log_target_density_new, log_q, log_target_density, log_q_new};

          // Rcout<<"rho="<<rho<<std::endl;
          // Rcout<<"test 8"<<std::endl;
          if(log(arma::randu())<=rho){
            theta_alpha(L_idx) = theta_alpha_new_block;
            alpha = alpha_new;
            accept(all_iter) = 1;
            accept_block(iter,m) = 1;
            temp = temp_new;
          }
          logll_M += -0.5*bma.square(arma::norm(temp,"fro"))*sigma_M2_inv;
  //         // Rcout<<"test 9"<<std::endl;
  //         // Rcout<<"idx="<<idx<<";idx%interval="<<idx%interval<< std::endl;
  //         // Rcout<<"iter%interval="<<iter%interval<<std::endl;
          if( (iter%interval==0) & (iter>0)  ){
            // Rcout<<"test 10"<<std::endl;
            arma::uvec u = arma::linspace<arma::uvec>(iter-interval,iter-1,interval);
            emp_accept.row(iter/interval-1) = mean(accept_block.rows(u),0);
            // Rcout<<"test1;mean(accept_block.rows(u),0)="<<mean(accept_block.rows(u),0)<<
            //   std::endl;
            // Rcout<<"test 11"<<std::endl;
            if(iter<stop_burnin){
              arma::colvec sigma_t = sqrt(2*step_all);
              for(arma::uword l = 0; l<num_block; l++){
                sigma_t(l) = bma.adjust_acceptance(emp_accept(iter/interval-1,l),sigma_t(l),target_accept);
                step_all(l) = sigma_t(l)*sigma_t(l)/2;
              }

            }

          }
  //         // when stop burnin, choose the average of last few steps
          if(iter==stop_burnin){
            // ivec back_temp = {interval*10,n_mcmc/2};
            // arma::uword back = min(back_temp);
            arma::uword back =  (n_mcmc - stop_burnin)/(n_mcmc/10) ;
            arma::uvec u = arma::linspace<arma::uvec>(iter-back,iter-1,back);
            step_all = mean(track_step.cols(u),1);
          }
        }// end of block update

        if( (iter%interval==0) & (iter>0)  ){
          // Rcout<<"test 10"<<std::endl;
          arma::uvec u = arma::linspace<arma::uvec>(iter-interval,iter-1,interval);
          emp_accept.row(iter/interval-1) = mean(accept_block.rows(u),0);
        }
      }// true end of block update
  //
  //     // Rcout<<"test 12"<<std::endl;
      track_step.col(iter) = step_all;
      all_iter = all_iter+1;
      track_grad_f.col(iter) = grad_f_all;
      track_rho_compo.slice(iter) = rho_compo;
      track_rho.col(iter) = rho_all;
      
      //     // -------------- Update lambda using RW MCMC --------------
      // double lambda_new = arma::randn()*sigma_lambda + lambda;
      // if(lambda_new>0 && lambda_new<1){
      //   double logll_M_new = 0;
      //   for(uword m=0; m<num_region; m++){
      //     arma::uvec delta = find(abs(alpha)>lambda_new);
      //     arma::uvec idx = region_idx[m];
      //     arma::uvec delta_Q = find(abs(alpha(idx))>lambda_new);
      //     mat Q = K[m];
      //     arma::uvec L_idx;
      //     if(m==0){
      //       L_idx = arma::linspace<uvec>(0,L_cumsum(m)-1,L_all(m));
      //     }else{
      //       L_idx = arma::linspace<uvec>(L_cumsum(m-1),L_cumsum(m)-1,L_all(m));
      //     }
      //     arma::uvec delta_in_block = intersect(idx,delta);
      // 
      //     mat K_block_t = Q.t();
      //     mat M_star = K_block_t*(M.rows(idx) - ones(idx.n_elem,1)*zetam.t()*C_t) - theta_eta.rows(L_idx);
      //     mat temp = M_star - K_block_t.cols(delta_Q)*(alpha.rows(delta_in_block)-sign(alpha.rows(delta_in_block))*lambda_new)*X.t();
      //     logll_M_new += -0.5*square(arma::norm(temp,"fro"))*sigma_M2_inv;
      //   }
      //   double rho = logll_M_new - logll_M;
      //   if(log(arma::randu()) <= rho){
      //     lambda = lambda_new;
      //     accept_lambda(iter) = 1;
      //   }
      // }
      // lambda_mcmc(iter) = lambda;
      // if( (iter%interval==0) & (iter>0) &(iter<stop_burnin) ){
      //   arma::uvec u = arma::linspace<arma::uvec>(iter-interval,iter-1,interval);
      //   double rate = arma::mean(accept.elem(u));
      //   sigma_lambda = adjust_acceptance(rate,sigma_lambda,0.4);
      // }
      

  //     // -------------- Update all other parameters using GS --------------
  arma::mat Mstar_alpha_term = arma::zeros(size(theta_eta));
          for(arma::uword m=0; m<num_region; m++){
            arma::uvec delta = find(abs(alpha)>lambda);
            arma::uvec idx = region_idx[m];
            arma::uvec delta_Q = find(abs(alpha(idx))>lambda);
            arma::uvec delta_in_block = intersect(idx,delta);
            arma::mat Q = K[m];
            arma::uvec L_idx;
            if(m==0){
              L_idx = arma::linspace< arma::uvec>(0,L_cumsum(m)-1,L_all(m));
            }else{
              L_idx = arma::linspace< arma::uvec>(L_cumsum(m-1),L_cumsum(m)-1,L_all(m));
            }
            arma::mat K_block_t = Q.t();
            Mstar_alpha_term.rows(L_idx) = K_block_t.cols(delta_Q)*(M.rows(delta_in_block)-( alpha.rows(delta_in_block)-sign(alpha.rows(delta_in_block))*lambda)*X.t() );
          }
          // 2. zetam
          arma::mat Sigma_zetam_inv = sigma_zetam2_inv + sigma_M2_inv*q2_sum*C_t*C;
          arma::mat Sigma_zetam = inv_sympd(Sigma_zetam_inv );
          arma::mat zeta_res = Mstar_alpha_term - theta_eta;
          arma::rowvec temp_zeta = q_vec.t()*zeta_res;
          arma::rowvec mu_zetam = temp_zeta*C *Sigma_zetam* sigma_M2_inv;
          zetam = arma::mvnrnd( mu_zetam.t(), Sigma_zetam);
          zetam_mcmc.col(iter) = zetam;
          
          
          // 3. theta_eta
          zeta_term =  q_vec*zetam.t()*C_t;
          arma::mat eta_res = Mstar_alpha_term - zeta_term; //L x n
          arma::colvec eta_Sigma_vec = 1/(sigma_M2_inv + sigma_eta2_inv/D);
          arma::colvec mu_eta_i;
          // mat mu_eta;
          for(arma::uword i=0;i<n;i++){
            mu_eta_i = (eta_res.col(i)*sigma_eta2_inv)%eta_Sigma_vec;
            arma::colvec theta_eta_i = arma::randn(L_max)%arma::sqrt(eta_Sigma_vec)+ mu_eta_i;
            theta_eta.col(i) = theta_eta_i-mean(theta_eta_i);
            // mu_eta = join_rows(mu_eta,mu_eta_i);
            // theta_eta.col(i) = mu_eta;
          }
          
          
          // 4. sigma_alpha
          sigma_alpha2_inv = arma::randg( arma::distr_param(a + L_max*0.5, 1/(b + dot(theta_alpha,theta_alpha/D)/2)) );
          sigma_alpha2_inv_mcmc(iter) = sigma_alpha2_inv;
          // 5. sigma_M
          arma::mat M_reg_res = Mstar_alpha_term - zeta_term - theta_eta;
          double M_norm = norm(M_reg_res,"fro");
          sigma_M2_inv = arma::randg( arma::distr_param(a + n*L_max/2,1/(b + M_norm*M_norm/2) ) );
          sigma_M2_inv_mcmc(iter) = sigma_M2_inv;
          
          //   6. sigma_zetam -> half-cauchy
          // sigma_zetam2_inv = arma::randg( arma::distr_param(a + zetam.n_elem/2,1/(b + dot(zetam,zetam)/2) ) );
          // sigma_zetam2_inv_mcmc(iter) = sigma_zetam2_inv;
          
          // //   7. sigma_eta -> giving too large variance
          double eta_norm = norm(theta_eta,"fro");
          sigma_eta2_inv = arma::randg( arma::distr_param(a + 0.5*n*L_max,1/(b + eta_norm*eta_norm/2)) );
          sigma_eta2_inv_mcmc(iter) = sigma_eta2_inv;
          // 
  
  //     // --------------------- summarize return --------------------- //
        
      theta_alpha_mcmc.col(iter) = theta_alpha;
      logll_mcmc(iter) = logll_M;
  }
  Rcpp::List gs = Rcpp::List::create(Rcpp::Named("theta_eta") = theta_eta,
                               Rcpp::Named("zetam_mcmc")= zetam_mcmc,
                               Rcpp::Named("lambda_mcmc")=lambda_mcmc,
                               Rcpp::Named("sigma_M2_inv_mcmc")= sigma_M2_inv_mcmc ,
                               Rcpp::Named("sigma_alpha2_inv_mcmc")= sigma_alpha2_inv_mcmc,
                               Rcpp::Named("sigma_zetam2_inv_mcmc")=  sigma_zetam2_inv_mcmc,
                               Rcpp::Named("sigma_eta2_inv_mcmc")= sigma_eta2_inv_mcmc);
  timer.step("end of iterations");
  return Rcpp::List::create(Rcpp::Named("theta_alpha_mcmc")=theta_alpha_mcmc,
                            Rcpp::Named("logll_mcmc")=logll_mcmc,
                            Rcpp::Named("track_grad_f")=track_grad_f,
                            Rcpp::Named("track_step")=track_step,
                            Rcpp::Named("accept_blcok")=accept_block,
                            Rcpp::Named("emp_accept")=emp_accept,
                            Rcpp::Named("track_rho")=track_rho,
                            Rcpp::Named("track_rho_compo")=track_rho_compo,
                            Rcpp::Named("gs")=gs,
                            Rcpp::Named("Timer")=timer
  );
}