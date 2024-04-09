# in development code for creating graphviz plots of GSEM models



fit_graph <- function(results, ...) {

  results_sort <- results %>% arrange(lhs, rhs)

  node_names <- unique(c(results_sort$lhs, results_sort$rhs))
  
  node_idx <- seq_along(node_names)
  names(node_idx) <- node_names
  
  edge_dir <- c("=~" = "forward", "~~" = "both")
  
  graph <- create_graph(
    nodes_df=create_node_df(n=length(node_names),
                            label=node_names, 
                            shape='oval', width=1,
                            #fillcolor=node_colors[node_names],
                            fontcolor='black'),
    edges_df=create_edge_df(from=node_idx[results_sort$lhs],
                            to=node_idx[results_sort$rhs],
                            label=round(results_sort$STD_Genotype, 2),
                            penwidth=0.3+abs(2*results_sort$STD_Genotype),
                            dir=edge_dir[results_sort$op]
                            ),
    attr_theme="tb")
  
  return(graph)

}

# render_fit <- function(results) render_graph(fit_graph(results))

# render_fit(pgc_commonfactor.fit$results)

# pgc_commonfactor.graph <- fit_graph(pgc_commonfactor.fit$results)

add_rank_same <- function(gv, top, bottom) {

  # add in block to specify node ranks
  gv.list <- str_split(gv, '\n\n')[[1]]
  
  # move the nodes/edges element to the end
  gv.list[6] <- gv.list[5]
  
  # add manually made ranks
  gv.list[5] <- paste("{rank=same",
  paste(paste0("'", top, "'"), collapse=' '),
  "}\n{rank=same",
  paste(paste0("'", bottom, "'"), collapse=' '),
  "}")
  
  rank.gv <- paste(gv.list, collapse='\n\n')

}


# pgc_commonfactor.gv <- add_rank_same(generate_dot(pgc_commonfactor.graph), 1, 2:5)
# grViz(pgc_commonfactor.gv)

# # output as a GraphViz dot file. Replace single quotes with double quotes 
# # as that's what the command line utility expects
# cat(str_replace_all(pgc_commonfactor_rank.gv, "'", '"'), file='mdd-symptom-gsem_files/pgc_commonfactor.gv')