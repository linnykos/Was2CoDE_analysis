
# per_cell_adjust = c("NB", "both"); quantile = 0.975; empirical_n = 1024

ideas_dist <-
  function(count_input, meta_cell, meta_ind, var_per_cell, var2test, 
           var2test_type = c("binary", "continuous"), 
           d_metric = c("Was", "JSD"), 
           fit_method = c("nb", "zinb", "kde", "dca_direct", "saver_direct"), 
           per_cell_adjust = c("NB", "both"), quantile = 0.975,
           empirical_n = 1024) {
    # -----------------------------------------------------------------
    # TZ:
    # initializes the var2test_type, d_metric, 
    # and fit_method parameters, validates the structure and content of 
    # meta_cell and meta_ind, which are 2 data frames,
    # checking for necessary columns like cell_id, individual, and var_per_cell
    # ----------------------------------------------------------------- 
    var2test_type   = var2test_type[1]
    d_metric        = d_metric[1]
    fit_method      = fit_method[1]
    
    
    if(!(is.data.frame(meta_cell))){
      stop("meta_cell should be a data.frame\n")
    }
    
    if(is.data.table(meta_cell)){
      meta_cell = as.data.frame(meta_cell)
    }
    
    if(!(is.data.frame(meta_ind))){
      stop("meta_ind should be a data.frame\n")
    }

    if(is.data.table(meta_ind)){
      meta_ind = as.data.frame(meta_ind)
    }

      # -----------------------------------------------------------------
      # TZ:
      # validates the structure and content of meta_cell and meta_ind, 
      # checking for necessary columns including cell_id, individual,
      # and the variable per cell (var_per_cell). 
      # ensures that these columns contain unique cell and individual IDs 
      # that match across the metadata and count matrices. 
      # aligning the gene expression data with the corresponding metadata 
      # for each cell and individual.
      # -----------------------------------------------------------------
      
      # -----------------------------------------------------------------
      # check the input data of count_input,when fit_method != dca_direct
      # -----------------------------------------------------------------
      count_matrix = count_input
      
      if(! is.matrix(count_matrix)){
        stop("count_matrix is not a matrix\n")
      }
      
      if( fit_method %in% c("nb", "zinb", "kde")){
        check_count <- function(v){any(v != round(v) | v < 0)}
        not_count = apply(count_matrix, 1, check_count)
        
        if(any(not_count)){
          str1 = "count_matrix should only include non-negative integers"
          str1 = sprintf("%s, violation in row %d\n", str1, which(not_count)[1])
          stop(str1)
        }
      }
      
      n_cell = ncol(count_matrix)
      n_gene = nrow(count_matrix)
      
      gene_ids = rownames(count_matrix)
      cell_ids = colnames(count_matrix)
      
      if(is.null(gene_ids)){
        stop("count_matrix should have row names for gene ids\n")
      }
      
      if(is.null(cell_ids)){
        stop("count_matrix should have col names for cell ids\n")
      }
      
      if(length(unique(gene_ids)) != n_gene){
        stop("row names of count_matrix (gene ids) are not unique\n")
      }
      
      if(length(unique(cell_ids)) != n_cell){
        stop("col names of count_matrix (cell ids) are not unique\n")
      }
      
      message(sprintf("the count_matrix includes %d genes in %d cells\n", 
                      n_gene, n_cell))
      
      # -----------------------------------------------------------------
      # check cell_id order of meta_cell, when fit_method != dca_direct
      # -----------------------------------------------------------------
      
      if(any(meta_cell$cell_id != colnames(count_matrix))){
        stop("cell_id in meta_cell do not match colnames of count_matrix\n")
      }
      
    

    # -----------------------------------------------------------------
    # check other aspects of meta_cell
    # -----------------------------------------------------------------    
    
    columns.meta.cell = c("cell_id", "individual", var_per_cell)
    
    if(! all(columns.meta.cell %in% names(meta_cell))){
      str1 = paste(columns.meta.cell, collapse=", ")
      stop(sprintf("names of meta_cell should contain %s\n", str1))
    }
    
    if(length(unique(meta_cell$cell_id)) != nrow(meta_cell)){
      stop("the cell_id's in meta_cell are not unique\n")
    }
    
    # -----------------------------------------------------------------
    # check the input data of meta_ind
    # -----------------------------------------------------------------
    
    columns.meta.ind = c("individual", var2test)
    
    if(! all(columns.meta.ind %in% names(meta_ind))){
      str1 = paste(columns.meta.ind, collapse=", ")
      stop(sprintf("names of meta_ind should contain %s\n", str1))
    }
    
    if(! setequal(meta_cell$individual, meta_ind$individual)){
      stop("the individual ids in meta_cell and meta_ind do not match\n")
    }
    
    if(length(unique(meta_ind$individual)) != nrow(meta_ind)){
      stop("the individual ids in meta_ind are not unique\n")
    }
    
    
    #if(any(meta_cell[,..var_per_cell] <= 0.0)){
    if(any(meta_cell[,var_per_cell, drop = FALSE] <= 0.0)){
      str1 = "the variables listed in 'var_per_cell' will be log transformed,"
      stop(paste(str1, "so they must be positive."))
    }
    
    
    # -----------------------------------------------------------------
    # TZ:
    # estimates the distribution(distance) for each gene and individual 
    # using Kernel Density Estimation (KDE)
    # -----------------------------------------------------------------
    
    # -----------------------------------------------------------------
    # estimate distance across individuals using kde
    # -----------------------------------------------------------------
    if (fit_method == "kde") {
      message("estimating distribution for each gene and each individual by kde\n")
      cov_value = apply(log10(meta_cell[,var_per_cell,drop=FALSE]), 2, median)
      #cov_value = apply(log10(meta_cell[, ..var_per_cell]), 2, median)
      dat_res=foreach (i_g = 1:n_gene) %dorng% {
        res_ig = list()
        # For loop 1:
        # For each gene (i_g), the function iterates over each individual 
        # in meta_ind, extracting the corresponding expression data 
        # from the count matrix. 
        # then applies a log transformation, 
        # preprocessing to make the data more normally distributed.
        # outputing the residuals that we want
        
        for (j in 1:nrow(meta_ind)) {
          ind_j = meta_ind$individual[j] #donor's ID
          w2use = which(meta_cell$individual == ind_j) 
          #grab all indexes associated with this donor
          base_j = c(t(lm_j$coefficients) %*% c(1, cov_value))
          res_ig[[j]] = lm_j$resid + base_j
        } #set up a linear reression log10(dat_j+const)~log10(library size)
        names(res_ig) = as.character(meta_ind$individual)
        res_ig #output the residual(*) and intercept
      }
      
      dist_array_list=foreach (i_g = 1:n_gene) %dorng% {
        
        res_ig = dat_res[[i_g]]
        dist_array1 = array(NA, dim=rep(nrow(meta_ind), 2))
        rownames(dist_array1) = meta_ind$individual
        colnames(dist_array1) = meta_ind$individual
        diag(dist_array1) = 0
        
        for (j_a in 1:(nrow(meta_ind)-1)) {
          res_a = res_ig[[j_a]]
          # For loop 2:
          # For each pair of donors 
          # compute the wasserstein distance b/w 2 distributions
                
          for (j_b in (j_a+1):nrow(meta_ind)) {
            res_b = res_ig[[j_b]]
            
            dist_array1[j_a, j_b] = tryCatch(
              divergence(res_a, res_b, d_metric = d_metric,  
                         fit_method = fit_method, empirical_n = empirical_n), #distance calculation 
              error = function(e) { NA }
            )
            
            dist_array1[j_b, j_a] = dist_array1[j_a, j_b]
          }
        }
        dist_array1
      }
    }
   
    # -----------------------------------------------------------------
    # TZ:
    # conclusion
    # compiles the pairwise distances between individuals for each gene 
    # into a 3-dimensional array (dist_array).  
    # -----------------------------------------------------------------
   
    length(dist_array_list)
    dim(dist_array_list[[1]])
    dist_array_list[[1]][1:2,1:2]
    
    nNA = sapply(dist_array_list, function(x){sum(is.na(c(x)))})
    table(nNA)
    
    dist_array = array(
      dim = c(
        n_gene,
        nrow(meta_ind),
        nrow(meta_ind)
      ),
      dimnames = list(gene_ids, meta_ind$individual, meta_ind$individual)
    )
    
    dim(dist_array)
    
    for (i in 1:n_gene){
      dist_array[i,,] = dist_array_list[[i]]
    }
    
    dim(dist_array)
    dist_array[1,1:2,1:2]
    
    dist_array
  }
