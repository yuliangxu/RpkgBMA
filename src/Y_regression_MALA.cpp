
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

#include "BMA_header.h"

// [[Rcpp::export]]
Rcpp::List Y_regression_GS(Rcpp::List& data, Rcpp::List& init,
                           Rcpp::List& region_idx_cpp,
                           Rcpp::List& kernel, arma::uword n_mcmc,
                           bool display_progress= true){
  
  Rcpp::Timer timer;  
  timer.step("start of precomputation"); 
  // input data
  arma::colvec Y = data["Y"];
  arma::colvec X = data["X"];
  arma::mat C = data["C"]; C = C.t();
  arma::mat M = data["M"];
  // input initial parameters
  arma::colvec zetay = init["zetay"];
  double sigma_beta = init["sigma_beta"]; 
  double sigma_Y = init["sigma_y"];
  double gamma = init["gamma"]; 
  double sigma_gamma = init["sigma_gamma"];
  double cy = init["cb"]; double a = init["a"]; double b = init["b"];
  double sigma_cy = init["sigma_cy"];
  double sigma_cy2 = sigma_cy*sigma_cy;
  double sigma_zeta_y = init["sigma_zeta_y"];
  double sigma_gamma2 = sigma_gamma*sigma_gamma;
  double sigma_beta2 = sigma_beta*sigma_beta;
  double sigma_Y2 = sigma_Y*sigma_Y;
  double sigma_zeta_y2 = sigma_zeta_y*sigma_zeta_y;
  
  // input basis functions
  Rcpp::List D = kernel["D"];Rcpp::List Q = kernel["Q"];
  
  
  // get M_star, D_vec, p, L from basis
  arma::uword num_region = region_idx_cpp.length();
  arma::uvec p_length; arma::uvec L_all;
  p_length.set_size(num_region); L_all.set_size(num_region);
  arma::mat M_star;
  arma::colvec D_vec;
  for(arma::uword i=0; i<num_region; i++){
    arma::mat Q_i = Q[i];
    p_length(i) = Q_i.n_rows; L_all(i)=Q_i.n_cols;
    arma::uvec  idx = region_idx_cpp[i]; 
    M_star = arma::join_cols(M_star, Q_i.t()*M.rows(idx));
    arma::colvec D_i = D[i];
    D_vec = arma::join_cols(D_vec, D_i);
  }
  arma::uword L = arma::sum(L_all); arma::uword p = arma::sum(p_length);
  M_star = M_star/sqrt(p);
  arma::mat M_star_t = M_star.t();
  
  arma::uword m = zetay.n_elem;
  arma::uvec L_cumsum = cumsum(L_all);
  
  // arma::colvec Y_star = Y- cy - gamma*X -  C*zetay;
  // arma::mat Sigma_theta_beta_inv = M_star* M_star_t/sigma_Y2 +
  //   diagmat(1/D_vec)/sigma_beta2;
  // arma::mat Sigma_theta_beta = inv_sympd(Sigma_theta_beta_inv);
  // arma::colvec mu_theta_beta = 1/sigma_Y2 * Sigma_theta_beta * M_star *
  //   Y_star;
  // colvec theta_beta = arma::mvnrnd( mu_theta_beta, Sigma_theta_beta );
  arma::colvec theta_beta = arma::ones(L,1);
  
  // initialize mcmc sequences
  arma::mat theta_beta_mcmc = arma::zeros(L,n_mcmc );
  arma::colvec gamma_mcmc = arma::zeros(n_mcmc,1);
  arma::mat zetay_mcmc = arma::zeros(zetay.n_elem,n_mcmc);
  arma::colvec cy_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_beta2_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_Y2_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_gamma2_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_zeta_y2_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_cy2_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec logLL_mcmc = arma::zeros(n_mcmc,1);
  
  Progress prog(n_mcmc*num_region, display_progress);
  timer.step("start of iteration"); 
  for(arma::uword iter=0; iter< n_mcmc; iter++){
    prog.increment(); 
    
    arma::colvec Y_star = Y- cy - gamma*X -  C*zetay;
    
    // 1. update theta_beta
    // for(uword m=0; m<num_region; m++){
    //   arma::uvec L_idx;
    //   if(m==0){
    //     L_idx = linspace<uvec>(0,L_cumsum(m)-1,L_all(m));
    //   }else{
    //     L_idx = linspace<uvec>(L_cumsum(m-1),L_cumsum(m)-1,L_all(m));
    //   }
    //   mat M_star_m = M_star.rows(L_idx);
    //   arma::mat Sigma_theta_beta_inv = M_star_m* M_star_m.t()/sigma_Y2 +
    //     diagmat(1/D_vec(L_idx))/sigma_beta2;
    //   arma::mat Sigma_theta_beta = inv_sympd(Sigma_theta_beta_inv);
    //   arma::colvec mu_theta_beta = 1/sigma_Y2 * Sigma_theta_beta * M_star_m *
    //     Y_star;
    //   theta_beta(L_idx) = arma::mvnrnd( mu_theta_beta, Sigma_theta_beta );
    // }
    arma::mat Sigma_theta_beta_inv = M_star* M_star_t/sigma_Y2 +
      diagmat(1/D_vec)/sigma_beta2;
    arma::mat Sigma_theta_beta = inv_sympd(Sigma_theta_beta_inv);
    arma::colvec mu_theta_beta = 1/sigma_Y2 * Sigma_theta_beta * M_star *
      Y_star;
    theta_beta = arma::mvnrnd( mu_theta_beta, Sigma_theta_beta );
    theta_beta_mcmc.col(iter) =  theta_beta;
    
    
    // 2. gamma
    double post_sigma_gamma2 = 1/(arma::sum(dot(X,X))/sigma_Y2 + 1/sigma_gamma2);
    double temp2 = dot(X, Y - cy - M_star_t*theta_beta - C * zetay);
    double mu_gamma = post_sigma_gamma2 * temp2/sigma_Y2;
    gamma = arma::randn() * sqrt(post_sigma_gamma2) + mu_gamma;
    gamma_mcmc(iter) = gamma;
    
    // 3. zeta_y
    arma::mat Sigma_zetay_inv = 1/sigma_Y2*(C.t()*C) + 1/sigma_zeta_y2;
    arma::mat Sigma_zetay = inv_sympd(Sigma_zetay_inv );
    arma::colvec mu_zetay = Sigma_zetay* C.t() *(  Y - cy - M_star_t*theta_beta - gamma*X )/sigma_Y2;
    zetay = arma::mvnrnd( mu_zetay, Sigma_zetay);
    zetay_mcmc.col(iter) = zetay;
    
    //   4. c_y 
    double post_sigma_cy = 1/(Y.n_elem/sigma_Y2 + 1/sigma_cy2);
    arma::colvec temp_cy = Y - M_star_t*theta_beta - C * zetay - gamma*X;
    double mu_cy = post_sigma_cy * arma::sum(temp_cy)/sigma_Y2;
    cy = arma::randn() * sqrt(post_sigma_cy) + mu_cy;
    cy_mcmc(iter) = cy;
    
    //   5. sigma_beta
    double sigma_beta_a = a + theta_beta.n_elem/2;
    double sigma_beta_b = b + dot(theta_beta,theta_beta/D_vec)/2;
    sigma_beta2 = 1/arma::randg( arma::distr_param(sigma_beta_a,1/sigma_beta_b) );
    sigma_beta2_mcmc(iter) = sigma_beta2;
    
    //   6. sigma_Y
    arma::colvec temp_sigma_Y = temp_cy-cy;
    double sigma_Y_b = b + dot(temp_sigma_Y,temp_sigma_Y)/2;
    sigma_Y2 = 1/arma::randg( arma::distr_param(a + Y.n_elem/2,1/sigma_Y_b) );
    sigma_Y2_mcmc(iter) = sigma_Y2;
    
    //   7. sigma_gamma -> use half-cauchy
    // sigma_gamma2 = 1/arma::randg( arma::distr_param(a + 0.5, 1/(b + gamma*gamma/2)) );
    // sigma_gamma2_mcmc(iter) = sigma_gamma2;
    
    //   8. sigma_zetay -> use half-cauchy
    // sigma_zeta_y2 = 1/arma::randg( arma::distr_param(a + zetay.n_elem/2,1/(b + dot(zetay,zetay)/2) ) );
    // sigma_zeta_y2_mcmc(iter) = sigma_zeta_y2;
    
    //   9. sigma_c_y -> use half-cauchy
    // sigma_cy2 = 1/arma::randg( arma::distr_param(a + 0.5,1/(b + cy*cy/2)) );
    // sigma_cy2_mcmc(iter) = sigma_cy2;
    
    // 10.logLL
    arma::colvec all_res = Y- cy - gamma*X -  C*zetay - M_star_t*theta_beta;
    logLL_mcmc(iter) = -0.5/sigma_Y2*dot(all_res,all_res);
  }
  timer.step("end of iterations");     
  
  return Rcpp::List::create(Rcpp::Named("theta_beta_mcmc")=theta_beta_mcmc,
                            Rcpp::Named("gamma_mcmc")=gamma_mcmc,
                            Rcpp::Named("zetay_mcmc")= zetay_mcmc,
                            Rcpp::Named("cy_mcmc")= cy_mcmc,
                            Rcpp::Named("timer")=timer,
                            Rcpp::Named("sigma_beta2_mcmc")= sigma_beta2_mcmc,
                            Rcpp::Named("sigma_Y2_mcmc")= sigma_Y2_mcmc ,
                            Rcpp::Named("sigma_gamma2_mcmc")= sigma_gamma2_mcmc,
                            Rcpp::Named("sigma_zeta_y2_mcmc")= sigma_zeta_y2_mcmc,
                            Rcpp::Named("sigma_cy2_mcmc")= sigma_cy2_mcmc,
                            Rcpp::Named("logLL_mcmc")= logLL_mcmc);
}

// -----------------------block update theta_beta--------------------------------- //
// [[Rcpp::export]]
Rcpp::List Y_regression_region_block_fast(arma::colvec& Y, arma::mat& M,
                                    arma::colvec& X, arma::mat& C,arma::uvec L_all,
                                    arma::uword num_region, Rcpp::List& region_idx,
                                    int n_mcmc, Rcpp::List& K, int stop_burnin,
                                    int unadjusted_langevin, int start_joint,
                                    double lambda, double target_accept,
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
  arma::uword L_max = arma::sum(L_all);
  // input
  arma::colvec theta_beta = init["theta_beta"];
  arma::colvec D = init["D"]; 
  arma::colvec D_sqrt = sqrt(D);
  double gamma = init["gamma"], cy = init["cb"]; 
  double a_sigma_beta = 1, b_sigma_beta = 1;
  arma::colvec zetay = init["zetay"];
  double sigma_Y = init["sigma_Y"], sigma_Y2 = sigma_Y*sigma_Y;
  double sigma_beta = init["sigma_beta"], sigma_beta2 = sigma_beta*sigma_beta;
  arma::colvec step_all = step*arma::ones(num_region);
  // hyper parameter for inverse-Gamma
  double a = 1.0, b=1.0, sigma_gamma2 = 1.0, mu_gamma = 1.0,sigma_zeta_y2 = 1.0,sigma_cy2=1.0;
  if(C.n_cols !=zetay.n_rows){
    Rcpp::Rcout<<"Error: dimensions of C and zetay don't match!"<<
      "dim of C = "<<size(C)<<"; dim of zetay = "<<size(zetay)<<std::endl;
    return Rcpp::List::create(Rcpp::Named("ERROR")=1);
  }
  arma::colvec Y_star = Y- cy - gamma*X -  C*zetay;
  arma::mat M_t = M.t(); arma::mat C_t = C.t();
  
  Rcpp::List Q_t(num_region) ;
  arma::colvec beta = arma::zeros(p,1);
  for(int l=0; l<num_region; l++){
    arma::uvec idx = region_idx[l];
    arma::mat Q = K[l];
    Q_t[l] = Q.t();
    arma::uvec L_idx;
    if(l==0){
      L_idx = arma::linspace<arma::uvec>(0,L_cumsum(l)-1,L_all(l));
    }else{
      L_idx = arma::linspace<arma::uvec>(L_cumsum(l-1),L_cumsum(l)-1,L_all(l));
    }
    beta(idx) = Q*theta_beta(L_idx);
  }
  
  //return 
  // arma::mat theta_beta_mcmc = zeros(L_max,n_mcmc);
  arma::mat theta_beta_mcmc_thin = arma::zeros(L_max,n_mcmc/interval);
  arma::colvec logll_mcmc_Y = arma::zeros(n_mcmc);
  // arma::mat track_grad_f = zeros(L_max,n_mcmc);
  // mat track_rho = zeros(num_region, n_mcmc);
  // mat track_step = zeros(num_region,n_mcmc);
  // cube track_rho_compo = zeros(num_region,4, n_mcmc);
  arma::mat emp_accept = arma::zeros( n_mcmc/interval,num_region);
  // arma::colvec accept = zeros(n_mcmc*num_region);
  arma::mat accept_block = arma::zeros(n_mcmc,num_region);
  // gs returns
  arma::colvec gamma_mcmc = arma::zeros(n_mcmc,1);
  arma::mat zetay_mcmc = arma::zeros(zetay.n_elem,n_mcmc);
  arma::colvec cy_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_beta2_mcmc = arma::zeros(n_mcmc,1);
  arma::colvec sigma_Y2_mcmc = arma::zeros(n_mcmc,1);
  // check if region_idx starts from 0!
  arma::uword min_region_idx=1;
  for(arma::uword l =0; l<num_region;l++){
    arma::uvec region_idx0 = region_idx[l];
    arma::uword min_l = min(region_idx0);
    if(min_l<min_region_idx){min_region_idx=min_l;}
  }
  
  if(min_region_idx>0){
    Rcpp::Rcout<<"Error: region_idx does not start from 0!"<<
      "min(region_idx[0]) = "<<min_region_idx<<std::endl;
    return Rcpp::List::create(Rcpp::Named("ERROR")=1,
                              Rcpp::Named("region_idx")=region_idx);
  }
  
  
  arma::uword all_iter=0;
  arma::uword num_block = num_region;
  arma::uvec delta_in_block;
  Progress prog(n_mcmc*num_region, display_progress);
  timer.step("start of iteration"); 
  for(int iter=0; iter<n_mcmc; iter++){
    // Rcout<<"iter="<<iter<<std::endl;
    if(iter==stop_burnin){
      // start the timer
      timer.step("stop of burnin");        // record the starting point
    }
    Y_star = Y- cy - gamma*X -  C*zetay;
    
    // colvec grad_f_all = zeros(L_max,1);
    // colvec rho_all = zeros(num_block,1);
    // mat rho_compo = zeros(num_block,4);
    
    
    // start block update
    double log_target_density;
    arma::colvec temp;
    for(arma::uword m=0; m < num_region; m++){
    prog.increment();
      arma::uvec delta = find(abs(beta)>lambda);
      arma::uvec idx = region_idx[m];
      arma::uvec delta_Q = find(abs(beta(idx))>lambda);
      arma::mat Q = K[m];
      arma::uvec L_idx;
      if(m==0){
        L_idx = arma::linspace<arma::uvec>(0,L_cumsum(m)-1,L_all(m));
      }else{
        L_idx = arma::linspace<arma::uvec>(L_cumsum(m-1),L_cumsum(m)-1,L_all(m));
      }

      // if(m==0){delta_in_block = delta_Q;}
      // else{ arma::uvec idx_ = region_idx[m-1];
      //   delta_in_block = delta_Q + idx_.n_elem;}
      if(m==0){
        delta_in_block = intersect(idx,delta);
      }

      arma::mat K_block_t = Q_t[m];
      arma::colvec temp_idx = M_t.cols(delta_in_block)*(beta.rows(delta_in_block)-sign(beta.rows(delta_in_block))*lambda)/sqrt(p);
      delta_in_block = intersect(idx,delta);



      if(m==0){
        temp = Y_star-M_t.cols(delta)*(beta.rows(delta)-sign(beta.rows(delta))*lambda)/sqrt(p);
      }

      arma::colvec grad_f = -theta_beta(L_idx)/D(L_idx)/sigma_beta2+
        K_block_t.cols(delta_Q)*M.rows(delta_in_block)*temp/sqrt(p)/sigma_Y2;
      double step = step_all(m);
      arma::colvec theta_beta_diff = step*grad_f+sqrt(2*step)*arma::randn(size(grad_f));
      arma::colvec theta_beta_new_block = theta_beta(L_idx)+theta_beta_diff;
      // grad_f_all(L_idx) = grad_f;



      if(iter< unadjusted_langevin){
        // accept(all_iter) = 1;
        accept_block(iter,m) = 1;
        beta(idx) = beta(idx) + Q*theta_beta_diff;
        theta_beta(L_idx) = theta_beta_new_block;

      }else{
        // MH step
        log_target_density = -0.5*bma.square(norm(theta_beta(L_idx)/D_sqrt(L_idx),2))/sigma_beta2-
          0.5*dot(temp,temp)/sigma_Y2;
        arma::colvec beta_new = beta;
        // beta_new(idx) = beta_new(idx) + Q*theta_beta_diff;
        beta_new(idx) += Q*theta_beta_diff;
        arma::uvec delta_new = find(abs(beta_new)>lambda);
        arma::uvec delta_Q_new = find(abs(beta_new(idx))>lambda);
        arma::uvec delta_in_block_new = intersect(idx,delta_new);

        // arma::colvec temp_new = Y_star-M_t.cols(delta_new)*(beta_new.rows(delta_new)-sign(beta_new.rows(delta_new))*lambda)/sqrt(p);
        arma::colvec x1 = arma::zeros(n,1); arma::colvec x2 = arma::zeros(n,1);
        bool b1 = delta_in_block.n_elem>0;
        bool b2 = delta_in_block_new.n_elem>0;
        if(b1){
          x1 = M_t.cols(delta_in_block)*(beta.rows(delta_in_block)-sign(beta.rows(delta_in_block))*lambda)/sqrt(p);
          // x1 = M_t.cols(delta_diff_in_block_old)*(beta.rows(delta_diff_in_block_old)-sign(beta.rows(delta_diff_in_block_old))*lambda)/sqrt(p);
        }
        if(b2){
          x2 = M_t.cols(delta_in_block_new)*(beta_new.rows(delta_in_block_new)-sign(beta_new.rows(delta_in_block_new))*lambda)/sqrt(p);
          // x2 = M_t.cols(delta_diff_in_block_new)*(beta_new.rows(delta_diff_in_block_new)-sign(beta_new.rows(delta_diff_in_block_new))*lambda)/sqrt(p);
            // - M_t.cols(delta_diff_in_block_new)*(beta.rows(delta_diff_in_block_new)-sign(beta.rows(delta_diff_in_block_new))*lambda)/sqrt(p);
        }
        arma::colvec temp_new;
        if(!b1 && b2){
          temp_new = temp-x2;
        }
        if(b1 && !b2){
          temp_new = temp+x1;
        }
        if(b1 && b2){
          temp_new = temp+x1-x2;
        }
        if(!b1 && !b2){
          temp_new = temp;
        }


        double log_target_density_new = -0.5*bma.square(norm( theta_beta_new_block/D_sqrt(L_idx),2))/sigma_beta2-
          0.5*dot(temp_new,temp_new)/sigma_Y2;
        arma::colvec grad_f_new = -theta_beta_new_block/D(L_idx)/sigma_beta2+
          K_block_t.cols(delta_Q_new)*M.rows(delta_in_block_new)*temp_new/sqrt(p)/sigma_Y2; // L x 1

        double log_q = -1/4/step * bma.square(norm(-theta_beta_diff-step*grad_f_new,2));
        double log_q_new = -1/4/step * bma.square(norm(theta_beta_diff-step*grad_f,2));
        double rho = log_target_density_new + log_q - log_target_density - log_q_new;
        // rho_all(m) = rho;
        // rho_compo.row(m) = {log_target_density_new, log_q, log_target_density, log_q_new};
        if(log(arma::randu())<=rho){
          theta_beta(L_idx) = theta_beta_new_block;
          beta = beta_new;
          temp = temp_new;
          log_target_density=log_target_density_new;
          // accept(all_iter) = 1;
          accept_block(iter,m) = 1;
        }
        if( (iter%interval==0) & (iter>0)  ){
          arma::uvec u = arma::linspace<arma::uvec>(iter-interval,iter-1,interval);
          emp_accept.row(iter/interval-1) = mean(accept_block.rows(u),0);
          if(iter<stop_burnin){
            arma::colvec sigma_t = sqrt(2*step_all);
            for(arma::uword l = 0; l<num_block; l++){
              sigma_t(l) = bma.adjust_acceptance(emp_accept(iter/interval-1,l),sigma_t(l),target_accept);
              step_all(l) = sigma_t(l)*sigma_t(l)/2;
            }

          }

        }
        // when stop burnin, choose the average of last few steps
        if(iter==stop_burnin){
          arma::uword back =  (n_mcmc - stop_burnin)/(n_mcmc/10) ;
          arma::uvec u = arma::linspace<arma::uvec>(iter-back,iter-1,back);
          
          // step_all = mean(track_step.cols(u),1);
        }
      }// end of block update
      if( (iter%interval==0) & (iter>0)  ){
        arma::uvec u = arma::linspace<arma::uvec>(iter-interval,iter-1,interval);
        emp_accept.row(iter/interval-1) = mean(accept_block.rows(u),0);
      }

      //----catch irregular updates
      // if(max(abs(theta_beta))>1e6 || abs(logll_mcmc_Y(iter))>1e6 ){
      //   Rcout<<"Error: max(abs(theta_beta))>1e6 || abs(logll_mcmc_Y(iter))>1e6 "<<std::endl;
      //   Rcpp::List result = Rcpp::List::create(Rcpp::Named("theta_beta_mcmc")=theta_beta_mcmc,
      //                                          // Rcpp::Named("accept")=accept,
      //                                          // Rcpp::Named("accept_rate")=emp_accept,
      //                                          Rcpp::Named("logll_mcmc_Y")=logll_mcmc_Y,
      //                                          Rcpp::Named("track_grad_f")=track_grad_f,
      //                                          Rcpp::Named("track_step")=track_step,
      //                                          Rcpp::Named("accept_blcok")=accept_block,
      //                                          Rcpp::Named("emp_accept")=emp_accept,
      //                                          Rcpp::Named("track_rho")=track_rho,
      //                                          Rcpp::Named("Timer")=timer);
      //   return Rcpp::List::create(Rcpp::Named("ERROR")=1,
      //                             Rcpp::Named("grad_f")=grad_f,
      //                             Rcpp::Named("step")=step,
      //                             Rcpp::Named("L_idx")=L_idx,
      //                             Rcpp::Named("idx")=idx,
      //                             Rcpp::Named("region")=m,
      //                             Rcpp::Named("delta_Q")=delta_Q,
      //                             Rcpp::Named("delta_in_block")=delta_in_block,
      //                             Rcpp::Named("temp")=temp,
      //                             Rcpp::Named("Y_star")=Y_star,
      //                             Rcpp::Named("delta")=delta,
      //                             Rcpp::Named("beta")=beta,
      //                             Rcpp::Named("result")=result);
      // }
    }


    // track_step.col(iter) = step_all;
    all_iter = all_iter+1;
    // track_grad_f.col(iter) = grad_f_all;
    // track_rho_compo.slice(iter) = rho_compo;
    // track_rho.col(iter) = rho_all;
    
    
    // -------------- Update GS for all other parameters --------------
    if(iter>start_joint){
      arma::uvec delta = find(abs(beta)>lambda);
      arma::colvec M_beta_term = M_t.cols(delta)*(beta.rows(delta)-sign(beta.rows(delta))*lambda)/sqrt(p);
      // // 2. gamma
      // double post_sigma_gamma2 = 1/(sum(dot(X,X))/sigma_Y2 + 1/sigma_gamma2);
      // double temp2 = dot(X, Y - cy - M_beta_term - C * zetay);
      // double mu_gamma = post_sigma_gamma2 * temp2/sigma_Y2;
      // gamma = arma::randn() * sqrt(post_sigma_gamma2) + mu_gamma;
      // gamma_mcmc(iter) = gamma;
      //
      // // 3. zeta_y
      arma::mat Sigma_zetay_inv = 1/sigma_Y2*(C_t*C) + 1/sigma_zeta_y2;
      arma::mat Sigma_zetay = inv_sympd(Sigma_zetay_inv );
      arma::colvec mu_zetay = Sigma_zetay* C_t *(  Y - cy - M_beta_term - gamma*X )/sigma_Y2;
      zetay = arma::mvnrnd( mu_zetay, Sigma_zetay);
      zetay_mcmc.col(iter) = zetay;
      //
      // //   4. c_y
      double post_sigma_cy = 1/(Y.n_elem/sigma_Y2 + 1/sigma_cy2);
      arma::colvec temp_cy = Y - M_beta_term - C * zetay - gamma*X;
      double mu_cy = post_sigma_cy * arma::sum(temp_cy)/sigma_Y2;
      cy = arma::randn() * sqrt(post_sigma_cy) + mu_cy;
      cy_mcmc(iter) = cy;
      //
      // //   5. sigma_beta
      double sigma_beta_a = a + theta_beta.n_elem/2;
      double sigma_beta_b = b + dot(theta_beta,theta_beta/D)/2;
      sigma_beta2 = 1/arma::randg( arma::distr_param(sigma_beta_a,1/sigma_beta_b) );
      sigma_beta2_mcmc(iter) = sigma_beta2;
      //
      // //   6. sigma_Y
      arma::colvec temp_sigma_Y = temp_cy-cy;
      double sigma_Y_b = b + dot(temp_sigma_Y,temp_sigma_Y)/2;
      sigma_Y2 = 1/arma::randg( arma::distr_param(a + Y.n_elem/2,1/sigma_Y_b) );
      sigma_Y2_mcmc(iter) = sigma_Y2;
      
      // 7. update lambda
      
    }

    // --------------------- summarize return --------------------- //
    
    // theta_beta_mcmc.col(iter) = theta_beta;
    arma::uvec delta_Y = find(abs(beta)>lambda);
    arma::colvec temp_Y = Y_star - M_t.cols(delta_Y)*(beta.rows(delta_Y)-sign(beta.rows(delta_Y))*lambda)/sqrt(p);
    logll_mcmc_Y(iter) = -dot(temp_Y,temp_Y)/2/sigma_Y2;
    
    if( (iter%interval==0) & (iter>0) ){
      theta_beta_mcmc_thin.col(iter/interval-1)=theta_beta;
    }
  }
  timer.step("end of iterations");   
  Rcpp::List gs = Rcpp::List::create(Rcpp::Named("gamma_mcmc")=gamma_mcmc,
                               Rcpp::Named("zetay_mcmc")= zetay_mcmc,
                               Rcpp::Named("cy_mcmc")= cy_mcmc,
                               Rcpp::Named("sigma_beta2_mcmc")= sigma_beta2_mcmc,
                               Rcpp::Named("sigma_Y2_mcmc")= sigma_Y2_mcmc);
  
  return Rcpp::List::create(
    // Rcpp::Named("theta_beta_mcmc")=theta_beta_mcmc,
                            Rcpp::Named("theta_beta_mcmc_thin")=theta_beta_mcmc_thin,
                            Rcpp::Named("logll_mcmc_Y")=logll_mcmc_Y,
                            // Rcpp::Named("track_grad_f")=track_grad_f,
                            // Rcpp::Named("track_step")=track_step,
                            Rcpp::Named("accept_blcok")=accept_block,
                            Rcpp::Named("emp_accept")=emp_accept,
                            // Rcpp::Named("track_rho")=track_rho,
                            // Rcpp::Named("track_rho_compo")=track_rho_compo,
                            Rcpp::Named("gs") = gs,
                            Rcpp::Named("Timer")=timer
  );
}