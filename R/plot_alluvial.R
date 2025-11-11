#' biowomp: Cluster-matching alluvial plots
#'
#' Main plotting function and helpers for bipartite-matching-based alluvial diagrams
#' @docType package
#' @name biowomp
#'
#' @importFrom dplyr mutate select group_by summarise arrange desc ungroup slice n pull filter bind_rows across matches all_of add_count distinct
#' @importFrom ggplot2 ggplot aes geom_text scale_fill_manual labs after_stat annotate theme_void theme element_text rel ggsave guides scale_color_manual scale_x_continuous element_blank coord_flip
#' @importFrom ggalluvial geom_alluvium geom_stratum stat_stratum stat_alluvium
#' @importFrom ggfittext geom_fit_text
#' @importFrom ggforce gather_set_data
#' @importFrom purrr map
#' @importFrom igraph max_bipartite_match V graph_from_data_frame cluster_louvain cluster_leiden E
#' @importFrom tibble is_tibble
#' @importFrom utils read.csv write.csv combn
#' @importFrom stats setNames
#' @importFrom rlang sym .data
#' @importFrom magrittr %>%
#' @importFrom data.table :=
#' @importFrom wompwomp data_preprocess data_sort
# devtools::load_all("~/Desktop/local/wompwomp")

utils::globalVariables(c(
    ".data", ":=", "%>%", "group_numeric", "col1_int", "col2_int", "id", "x", "y", "value", "stratum", "total", "cum_y", "best_cluster_agreement", "neighbor_net", "alluvium", "pos", "count", "group1", "group2", "value", "group1_size", "group2_size", "weight", "parent", "group_name"
))

StatStratum <- ggalluvial::StatStratum # avoid the error Can't find stat called "stratum" - and make sure to do stat = StatStratum instead of stat = "stratum"

default_colors <- c(
    "#D55E00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#E69F00", "#CC79A7", "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711", "#005685",
    "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF", "#00F6B3", "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C", "#FFCB57", "#9AD2F2",
    "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D", "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D", "#00446B", "#803800", "#8D3666",
    "#3D3D3D"
)

compute_alluvial_statistics <- function(clus_df_gather, graphing_columns, column_weights = "value") {
    message(sprintf("Alluvial statistics: n = number of elements; m = number of graphing columns; a = number of alluvia/edges; k_i = number of blocks in layer i (where i goes from 1:m); K_sum = number of blocks across all layers; K_prod = product of blocks across all layers"))
    message(sprintf("n = %s", sum(clus_df_gather[[column_weights]], na.rm = TRUE)))
    message(sprintf("m = %s", length(graphing_columns)))
    message(sprintf("a = %s", nrow(clus_df_gather)))
    
    K_sum <- 0
    K_prod <- 1
    for (i in 1:length(graphing_columns)) {
        colname <- paste0("col", i, "_int")
        k_i <- length(unique(clus_df_gather[[colname]]))
        message(sprintf("k_%s = %s", i, k_i))
        
        K_sum <- K_sum + k_i
        K_prod <- K_prod * k_i
    }
    
    message(sprintf("K_sum = %s", K_sum))
    message(sprintf("K_prod = %s", K_prod))
}

load_in_df <- function(df, graphing_columns = NULL, column_weights = NULL) {
    if (is.character(df) && grepl("\\.csv$", df)) {
        df <- read.csv(df) # load in CSV as dataframe
    } else if (tibble::is_tibble(df)) {
        df <- as.data.frame(df) # convert tibble to dataframe
    } else if (!is.data.frame(df)) {
        stop("Input must be a data frame, tibble, or CSV file path.")
    }
    
    if (!is.null(column_weights)) {
        if (!(column_weights %in% colnames(df))) {
            stop(sprintf("column_weights '%s' is not a column in the dataframe.", column_weights))
        }
        df <- tidyr::uncount(df, weights = !!rlang::sym(column_weights))
    }
    
    if (!(is.null(graphing_columns))) {
        for (col in graphing_columns) {
            if (!(col %in% colnames(df))) {
                stop(sprintf("column '%s' is not a column in the dataframe.", col))
            }
            # convert to factor
            if (!is.factor(df[[col]])) {
                df[[col]] <- as.factor(df[[col]])
            }
        }
    }
    
    return(df)
}

lowercase_args <- function(arg_names) {
    for (nm in arg_names) {
        val <- get(nm, envir = parent.frame())
        if (is.character(val)) {
            assign(nm, tolower(val), envir = parent.frame())
        }
    }
}

print_function_params <- function() {
    f <- sys.function(sys.parent())
    call <- sys.call(sys.parent())
    defaults <- as.list(formals(f))
    call_args <- as.list(call)[-1]
    
    # Evaluate args, but wrap in list() so NULL survives
    eval_args <- lapply(call_args, function(arg) {
        if (is.symbol(arg) && as.character(arg) == "NULL") {
            list(NULL)  # wrapper preserves NULL
        } else {
            list(eval(arg, envir = parent.frame()))
        }
    })
    
    # Flatten one level
    names(eval_args) <- names(call_args)
    eval_args <- lapply(eval_args, `[[`, 1)
    
    # Manual merge (defaults first, then override)
    all_args <- defaults
    for (nm in names(eval_args)) {
        all_args[nm] <- list(eval_args[[nm]])  # assign inside list()
    }
    
    # Print
    for (nm in names(all_args)) {
        val <- all_args[[nm]]
        if (is.null(val)) {
            message(nm, " = NULL")
        } else if (length(val) == 1) {
            message(nm, " = ", val)
        } else {
            message(nm, " =")
            for (j in seq_along(val)) {
                elt_name <- names(val)[j]
                if (!is.null(elt_name) && nzchar(elt_name)) {
                    message("  ", elt_name, " : ", val[[j]])
                } else {
                    message("  [", j, "] : ", val[[j]])
                }
            }
        }
    }
}

find_colors_advanced <- function(clus_df_gather, ditto_colors, graphing_columns, coloring_algorithm_advanced_option = "leiden", resolution = 1) {
    column_int_names <- c()
    for (col_int in seq_along(graphing_columns)) {
        int_name <- paste0("col", col_int, "_int")
        column_int_names <- c(column_int_names, int_name)
    }
    clus_df_ungrouped <- clus_df_gather[, c(column_int_names, "value")]
    
    first <- TRUE
    compared <- c()
    for (group1_name in column_int_names) {
        for (group2_name in column_int_names) {
            if (!(group1_name == group2_name)) {
                comp1 <- paste0(group1_name, group2_name)
                comp2 <- paste0(group2_name, group1_name)
                if (!(comp1 %in% compared | comp2 %in% compared)) {
                    if (first) {
                        clus_df_filtered <- clus_df_ungrouped[, c(group1_name, group2_name, "value")]
                        clus_df_filtered <- clus_df_filtered %>%
                            add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) %>%
                            select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
                        clus_df_filtered <- distinct(clus_df_filtered)
                        colnames(clus_df_filtered) <- c("group1", "group2", "value")
                        
                        clus_df_filtered <- clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(group1_size = sum(value))
                        clus_df_filtered <- clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(group2_size = sum(value))
                        
                        clus_df_filtered <- clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(weight = value) # /group2_size)
                        
                        clus_df_filtered$group1 <- sub("^", paste0(group1_name, "_"), clus_df_filtered[["group1"]])
                        clus_df_filtered$group2 <- sub("^", paste0(group2_name, "_"), clus_df_filtered[["group2"]])
                        
                        first <- FALSE
                    } else {
                        temp_clus_df_filtered <- clus_df_ungrouped[, c(group1_name, group2_name, "value")]
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) %>%
                            select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
                        temp_clus_df_filtered <- distinct(temp_clus_df_filtered)
                        colnames(temp_clus_df_filtered) <- c("group1", "group2", "value")
                        
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(group1_size = sum(value))
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            group_by(group2) %>%
                            mutate(group2_size = sum(value))
                        
                        temp_clus_df_filtered <- temp_clus_df_filtered %>%
                            group_by(group1) %>%
                            mutate(weight = value) # /group2_size)
                        
                        temp_clus_df_filtered$group1 <- sub("^", paste0(group1_name, "_"), temp_clus_df_filtered[["group1"]])
                        temp_clus_df_filtered$group2 <- sub("^", paste0(group2_name, "_"), temp_clus_df_filtered[["group2"]])
                        
                        clus_df_filtered <- rbind(clus_df_filtered, temp_clus_df_filtered)
                    }
                }
                # compared <- c(compared, comp1, comp2)
            }
        }
    }
    
    clus_df_extra_filtered <- clus_df_filtered[, c("group1", "group2", "value")]
    g <- igraph::graph_from_data_frame(d = clus_df_extra_filtered, directed = FALSE)
    if (coloring_algorithm_advanced_option == "louvain") {
        partition <- igraph::cluster_louvain(g, weights = igraph::E(g)$value, resolution = resolution)
    } else if (coloring_algorithm_advanced_option == "leiden") {
        partition <- igraph::cluster_leiden(g, weights = igraph::E(g)$value, resolution = resolution)
    } else {
        stop(sprintf("coloring_algorithm_advanced_option '%s' is not recognized. Please choose from 'leiden' (default) or 'louvain'.", coloring_algorithm_advanced_option))
    }
    
    clus_df_leiden <- data.frame(group_name = partition$names, leiden = partition$membership)
    clus_df_leiden <- clus_df_leiden %>% tidyr::separate_wider_delim(group_name, names = c("col", "trash", "group"), delim = "_")
    clus_df_leiden <- clus_df_leiden %>% tidyr::separate_wider_delim(col, names = c("trash2", "col"), delim = "col")
    clus_df_leiden <- clus_df_leiden[, c("col", "group", "leiden")]
    clus_df_leiden[["colors"]] <- unlist(Map(function(x) ditto_colors[x], clus_df_leiden$leiden))
    
    final_colors <- c()
    for (col_int in unique(clus_df_leiden$col)) {
        temp_df <- clus_df_leiden[clus_df_leiden$col == col_int, ]
        final_colors <- c(final_colors, temp_df[rev(order(temp_df$group)), ]$colors)
    }
    
    return(final_colors)
}

find_group2_colors <- function(clus_df_gather, ditto_colors, unused_colors, current_g1_colors,
                               group1_name = "col1_int", group2_name = "col2_int",
                               cutoff = .5) {
    num_levels <- length(levels(clus_df_gather[[group2_name]]))
    
    clus_df_ungrouped <- clus_df_gather[, c(group1_name, group2_name, "value")]
    clus_df_filtered <- clus_df_ungrouped %>%
        add_count(!!rlang::sym(group1_name), !!rlang::sym(group2_name), wt = value) %>%
        select(!!rlang::sym(group1_name), !!rlang::sym(group2_name), n)
    clus_df_filtered <- distinct(clus_df_filtered)
    colnames(clus_df_filtered) <- c(group1_name, group2_name, "value")
    
    
    clus_df_filtered <- clus_df_filtered %>%
        group_by(!!rlang::sym(group1_name)) %>%
        mutate(group1_size = sum(value))
    clus_df_filtered <- clus_df_filtered %>%
        group_by(!!rlang::sym(group2_name)) %>%
        mutate(group2_size = sum(value))
    clus_df_filtered <- clus_df_filtered %>%
        group_by(!!rlang::sym(group1_name)) %>%
        mutate(weight = value / group2_size)
    
    parent_df <- clus_df_filtered %>%
        group_by(!!rlang::sym(group2_name)) %>%
        filter(weight == max(weight)) # %>% select(!!rlang::sym(group1_name))
    colnames(parent_df)[which(names(parent_df) == group1_name)] <- "parent"
    parent_df <- parent_df %>% mutate(parent = (ifelse(weight > cutoff, as.character(parent), "")))
    parent_df <- parent_df[, c(group2_name, "parent")]
    clus_df_filtered <- merge(clus_df_filtered, parent_df, by = group2_name)
    
    clus_df_filtered[["colors"]] <- unlist(Map(function(x) ifelse(x == "", "", current_g1_colors[x]), clus_df_filtered$parent))
    clus_df_parent <- clus_df_filtered[clus_df_filtered$parent != "", ]
    group2_colors <- vector("character", num_levels)
    for (col_int in unique(clus_df_gather[[group2_name]])) {
        if (col_int %in% unique(clus_df_parent[[group2_name]])) {
            temp_df <- clus_df_parent[clus_df_parent[[group2_name]] == col_int, ]
            temp_df <- temp_df[temp_df[[group1_name]] == temp_df[["parent"]], ]
            group2_colors[[as.integer(col_int)]] <- temp_df[["colors"]]
        } else {
            group2_colors[[as.integer(col_int)]] <- unused_colors[1]
            unused_colors <- unused_colors[2:length(unused_colors)]
        }
    }
    
    return(group2_colors)
}


#' Generate an Alluvial Plot with Minimal Cluster Cross-over
#'
#' Creates a two-axis alluvial plot to visualize the relationship between two categorical groupings (e.g., cluster assignments from different methods),
#' minimizing crossover and aligning clusters based on agreement.
#'
#' @param df A data frame, tibble, or CSV file path. Must be in one of two formats:
#' (1) column_weights == NULL: Each row represents an entity, each column represents a grouping, and each entry represents the membership of the entity in that row to the grouping in that column. Must contain at least two columns (two graphing_columns).
#' (2) column_weights != NULL: Each row represents a combination of groupings, each column from \code{graphing_columns} represents a grouping, and the column \code{column_weights} represents the number of entities in that combination of groupings. Must contain at least three columns (two \code{graphing_columns}, one \code{column_weights}).
#' @param graphing_columns Optional character vector. Vector of column names from \code{df} to be used in graphing (i.e., alluvial plotting). Mutually exclusive with \code{column1} and \code{column2}.
#' @param column1 Optional character. Can be used along with \code{column2} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column2 Optional character. Can be used along with \code{column1} in place of \code{graphing_columns} if working with two columns only. Mutually exclusive with \code{graphing_columns}.
#' @param column_weights Optional character. Column name from \code{df} that contains the weights of each combination of groupings if \code{df} is in format (2) (see above).
#' @param sorting_algorithm Character. Algorithm with which to sort the values in the dataframe. Can choose from: 'neighbornet', 'tsp', 'greedy_wolf', 'greedy_wblf', 'random', 'none'. 'neighbornet' performs sorting with NeighborNet (Bryant and Moulton, 2004). 'tsp' performs Traveling Salesman Problem solver from the TSP package. 'greedy_wolf' implements a custom greedy algorithm where one layer is fixed, and the other layer is sorted such that each node is positioned as close to its largest parent from the fixed side as possible in a greedy fashion. 'greedy_wblf' implements the 'greedy_wolf' algorithm described previously twice, treating each column as fixed in one iteration and free in the other iteration. 'greedy_wolf' and 'greedy_wblf' are only valid when \code{graphing_columns} has exactly two entries. 'random' randomly maps blocks. 'none' keeps the mappings as-is when passed into the function.
#' @param optimize_column_order Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param optimize_column_order_per_cycle Logical. If TRUE, will optimize the order of \code{graphing_columns} to minimize edge overlap upon each cycle. If FALSE, will optimize the order of \code{graphing_columns} to minimize edge overlap on the beginning cycle only. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{length(graphing_columns) > 2}.
#' @param matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in different layers without a shared edge/path. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param same_side_matrix_initialization_value Positive integer. Initialized value in distance matrix for nodes in the same layer. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param weight_scalar Positive integer. Scalar with which to multiply edge weights after taking their -log in the distance matrix for nodes with a nonzero edge. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param matrix_initialization_value_column_order Positive integer. Initialized value in distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param weight_scalar_column_order Positive integer. Scalar with which to loss function after taking their log1p in the distance matrix for optimizing column order. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_metric Character. Metric to use for determining column order. Options are "edge_crossing" (default) or "ARI". Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param column_sorting_algorithm Character. Algorithm to use for determining column order. Options are "tsp" (default) or "neighbornet". Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'} and \code{optimize_column_order} is TRUE.
#' @param weighted Logical. Weighted objective
#' @param cycle_start_positions Set. Cycle start positions to consider. Anything outside this set will be skipped. Only applies when \code{sorting_algorithm == 'neighbornet' or 'tsp'}.
#' @param fixed_column Character or Integer. Name or position of the column in \code{graphing_columns} to keep fixed during sorting. Only applies when \code{sorting_algorithm == 'greedy_wolf'}.
#' @param random_initializations Integer. Number of random initializations for the positions of each grouping in \code{graphing_columns}. Only applies when \code{sorting_algorithm == 'greedy_wolf' or sorting_algorithm == 'greedy_wblf'}.
#' @param color_boxes Logical. Whether to color the strata/boxes (representing groups).
#' @param color_bands Logical. Whether to color the alluvia/edges (connecting the strata).
#' @param color_list Optional named list or vector of colors to override default group colors.
#' @param color_band_list Optional named list or vector of colors to override default band colors.
#' @param color_band_column Optional Character. Which column to use for coloring bands.
#' @param color_val Optional named list where the entries are colors and the names correspond to values of the dataframe that should use those colors
#' @param color_band_boundary Logical. Whether or not to color boundaries between bands
#' @param coloring_algorithm Character. Matching colors methods. Choices are 'advanced' (default), 'none', 'left', 'right', or any value in \code{graphing_columns}.
#' @param coloring_algorithm_advanced_option Character. If \code{coloring_algorithm == 'advanced'}, then choose graph clustering algorithm. Choices are 'leiden' (default) or 'louvain'.
#' @param resolution Numeric If \code{coloring_algorithm == 'advanced'}, then choose resolution for the graph clustering algorithm. Affects coloring of both bands and boxes.
#' @param cutoff Numeric If \code{coloring_algorithm != 'none' and coloring_algorithm != 'advanced'}, sets the cutoff for color matching, below which a new color will be assigned.
#' @param alluvial_alpha Numeric between 0 and 1. Transparency level for the alluvial bands.
#' @param include_labels_in_boxes Logical. Whether to include text labels inside the rectangular group boxes.
#' @param include_axis_titles Logical. Whether to display axis titles for column1 and column2.
#' @param include_group_sizes Logical. If \code{TRUE}, includes group sizes in the labels (e.g., "Group A (42)").
#' @param output_plot_path Character. File path to save the plot (e.g., "plot.png"). If \code{NULL}, then will not be saved.
#' @param output_df_path Optional character. Output path for the output data frame, in CSV format. If \code{NULL}, then will not be saved.
#' @param preprocess_data Logical. If TRUE, will preprocess the data with the \code{data_preprocess} function.
#' @param default_sorting Character. Default column sorting in data_preprocess if integer columns do not exist. Options are 'alphabetical' (default), 'reverse_alphabetical', 'increasing', 'decreasing', 'random'.
#' @param box_width Numeric between 0 and 1. Box width
#' @param text_width Numeric between 0 and 1. Text width
#' @param min_text Integer greater than 0. Min text
#' @param text_size Integer greater than 0. Text size (works whether auto_adjust_text is TRUE or FALSE).
#' @param auto_adjust_text Logical. Whether to automatically adjust text size to fit in box.
#' @param axis_text_size Integer greater than 0. Axis text size
#' @param axis_text_vjust Integer. Axis text vjust
#' @param save_height Integer greater than 0. Save height, in inches
#' @param save_width Integer greater than 0. Save width, in inches
#' @param dpi Integer greater than 0. DPI for \code{output_plot_path}, if \code{output_plot_path} is a raster image or \code{rasterise_alluvia} is TRUE
#' @param rasterise_alluvia Logical. Whether to rasterize the alluvia if \code{output_plot_path} is a PDF. Can save space if DPI low enough
#' @param keep_y_labels Keep y labels
#' @param keep_x_labels Keep x labels
#' @param box_line_width Box line width
#' @param verbose Logical. If TRUE, will display messages during the function.
#' @param print_params Logical. If TRUE, will print function params.
#' @param add_legend Logical. If TRUE, will generate a legend of the colors of boxes and alluvial
#' @param legend_loc Character. Location of legend. Only applies if \code{add_legened == TRUE}. Choices are 'right' (default), 'left', 'bottom', 'top'
#' @param flip_xy Logical. Flip x and y (rotate plot 90 degrees).
#'
#' @return A \code{ggplot2} object representing the alluvial plot.
#'
#' @examples
#' # Example 1: df format 1
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' p <- plot_alluvial(df,
#'     graphing_columns = c("method1", "method2"),
#'     sorting_algorithm = "tsp",
#'     column_sorting_algorithm = "tsp",
#'     coloring_algorithm = "right"
#' )
#'
#' # Example 2: df format 2
#' df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
#' clus_df_gather <- df |>
#'     dplyr::mutate_if(is.numeric, function(x) factor(x, levels = as.character(sort(unique(x))))) |>
#'     dplyr::group_by_all() |>
#'     dplyr::count(name = "value")
#' p <- plot_alluvial(
#'     clus_df_gather,
#'     graphing_columns = c("method1", "method2"),
#'     column_weights = "value",
#'     sorting_algorithm = "tsp",
#'     column_sorting_algorithm = "tsp",
#'     coloring_algorithm = "right"
#' )
#'
#' @export
plot_alluvial <- function(df, graphing_columns = NULL, column1 = NULL, column2 = NULL,
                          column_weights = NULL, sorting_algorithm = "tsp",
                          optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE,
                          matrix_initialization_value = 1e6, same_side_matrix_initialization_value = 1e6,
                          weight_scalar = 5e5, matrix_initialization_value_column_order = 1e6,
                          weight_scalar_column_order = 1, column_sorting_metric = "edge_crossing",
                          column_sorting_algorithm = "tsp", weighted = TRUE, cycle_start_positions = NULL, fixed_column = NULL,
                          random_initializations = 1, color_boxes = TRUE, color_bands = FALSE,
                          color_list = NULL, color_band_list = NULL, color_band_column = NULL, color_val = NULL,
                          color_band_boundary = FALSE, coloring_algorithm = "advanced", coloring_algorithm_advanced_option = "leiden", resolution = 1, cutoff = .5,
                          alluvial_alpha = 0.5, include_labels_in_boxes = TRUE, include_axis_titles = TRUE, include_group_sizes = FALSE,
                          output_plot_path = NULL, output_df_path = NULL, preprocess_data = TRUE,
                          default_sorting = "alphabetical", box_width = 1 / 3, text_width = 1 / 4, min_text = 4, text_size = 14,
                          auto_adjust_text = TRUE, axis_text_size = 2, axis_text_vjust = 0, save_height = 6, save_width = 6, dpi = 300, rasterise_alluvia = FALSE,
                          keep_y_labels = FALSE, keep_x_labels = TRUE, 
                          box_line_width = 1, verbose = FALSE, print_params = FALSE,
                          add_legend = FALSE, legend_loc = "right", flip_xy=FALSE) {
    if (print_params) print_function_params()
    lowercase_args(c("sorting_algorithm", "column_sorting_metric", "column_sorting_algorithm", "coloring_algorithm", "coloring_algorithm_advanced_option", "default_sorting", "legend_loc"))
    
    #* Type Checking Start
    # ensure someone doesn't specify both graphing_columns and column1/2
    if (!is.null(graphing_columns) && (!is.null(column1) || !is.null(column2))) {
        stop("Specify either graphing_columns or column1/column2, not both.")
    }
    
    if (is.vector(graphing_columns) && length(graphing_columns) < 2) {
        stop("graphing_columns must have at least 2 entries.")
    }
    
    if (preprocess_data) {
        if (verbose) message("Loading in data")
        df <- load_in_df(df = df, graphing_columns = graphing_columns, column_weights = column_weights)
    }
    
    if (nrow(df) == 0) {
        stop("df has no rows.")
    }
    
    if (!is.null(graphing_columns) && any(!graphing_columns %in% colnames(df))) {
        stop("Some graphing_columns are not present in the dataframe.")
    }
    
    if (ncol(df) < 2) {
        stop("Dataframe must have at least 2 columns when column_weights is NULL.")
    } else if (ncol(df) > 2) {
        if (is.null(graphing_columns) && is.null(column1) && is.null(column2)) {
            stop("graphing_columns must be specified when dataframe has more than 2 columns and column_weights is NULL.")
        }
    } else { # length 2
        if (is.null(column1) && !is.null(column2)) {
            column1 <- setdiff(colnames(df), column2)
        } else if (is.null(column2) && !is.null(column1)) {
            column2 <- setdiff(colnames(df), column1)
        } else if (is.null(column1) && is.null(column2)) {
            column1 <- colnames(df)[1]
            column2 <- colnames(df)[2]
        }
    }
    
    # if someone specifies column1/2, then use it
    if (length(graphing_columns) == 2) {
        column1 <- graphing_columns[1]
        column2 <- graphing_columns[2]
    }
    
    if (is.null(graphing_columns)) {
        graphing_columns <- c(column1, column2)
    }
    
    if (is.null(fixed_column)) {
        fixed_column <- column1
    } else if ((is.integer(fixed_column) || (is.double(fixed_column)))) {
        if (fixed_column > length(colnames(df))) {
            stop(sprintf("fixed_column index '%s' is not a column in the dataframe.", fixed_column))
        } else {
            fixed_column <- colnames(df)[fixed_column]
        }
    } else if (!(fixed_column %in% colnames(df))) {
        stop(sprintf("fixed_column '%s' is not a column in the dataframe.", fixed_column))
    }
    
    if (!is.null(color_band_column)) {
        color_bands <- TRUE
    }
    #* Type Checking End
    
    if (sorting_algorithm == "greedy_wolf") {
        default_sorting <- "fixed"
    }
    # Preprocess
    if (preprocess_data) {
        if (verbose) message("Preprocessing data before sorting")
        clus_df_gather_unsorted <- wompwomp::data_preprocess(df = df, graphing_columns = graphing_columns, column_weights = column_weights, color_band_column = color_band_column, default_sorting = default_sorting, load_df = FALSE, do_gather_set_data = FALSE)
        if (is.null(column_weights)) {
            column_weights <- "value" # is set during data_preprocess
        }
    } else {
        clus_df_gather_unsorted <- df
    }
    
    do_compute_alluvial_statistics <- TRUE
    if (verbose) {
        compute_alluvial_statistics(clus_df_gather = clus_df_gather_unsorted, graphing_columns = graphing_columns, column_weights = column_weights)
        do_compute_alluvial_statistics <- FALSE
    }
    
    # Sort
    if (verbose) message(sprintf("Sorting data with sorting_algorithm=%s", sorting_algorithm))
    data_sort_output <- wompwomp::data_sort(df = clus_df_gather_unsorted, graphing_columns = graphing_columns, column_weights = column_weights, sorting_algorithm = sorting_algorithm, optimize_column_order = optimize_column_order, optimize_column_order_per_cycle = optimize_column_order_per_cycle, matrix_initialization_value = matrix_initialization_value, same_side_matrix_initialization_value = same_side_matrix_initialization_value, weight_scalar = weight_scalar, weighted = weighted, matrix_initialization_value_column_order = matrix_initialization_value_column_order, weight_scalar_column_order = weight_scalar_column_order, column_sorting_metric = column_sorting_metric, column_sorting_algorithm = column_sorting_algorithm, cycle_start_positions = cycle_start_positions, fixed_column = fixed_column, random_initializations = random_initializations, output_df_path = output_df_path, return_updated_graphing_columns = TRUE, preprocess_data = FALSE, load_df = FALSE, verbose = verbose, do_compute_alluvial_statistics = do_compute_alluvial_statistics)
    df <- data_sort_output$clus_df_gather  #!!! I'm overriding df here
    graphing_columns <- data_sort_output$graphing_columns
    
    # Plot
    if (verbose) message("Plotting data")
    #* beginning of the old plot_alluvial_internal function
    lowercase_args(c("coloring_algorithm", "coloring_algorithm_advanced_option", "legend_loc"))
    
    if (verbose) compute_alluvial_statistics(clus_df_gather = df, graphing_columns = graphing_columns, column_weights = column_weights)
    
    geom_alluvium <- if (rasterise_alluvia) {
        function(...) ggrastr::rasterise(ggalluvial::geom_alluvium(...), dpi = dpi)
    } else {
        ggalluvial::geom_alluvium
    }
    
    if (column_weights != "value") {
        df <- df %>% dplyr::rename(value = !!sym(column_weights))
        column_weights <- "value"
    }
    
    df <- ggforce::gather_set_data(df, 1:2)
    if (!is.numeric(df$x)) {
        df$x <- match(as.character(df$x), graphing_columns)
    } # weird Docker issue
    df <- df[df$x == 1, ]
    
    if (!is.null(color_list)) {
        ditto_colors <- color_list
    } else {
        ditto_colors <- default_colors
    }
    
    # remove user-defined colors from available color list
    if (!(is.null(color_val))) {
        # convert named list into named vector
        if (is.list(color_val)) {
            color_val <- unlist(color_val)
        }
        ditto_colors <- ditto_colors[!(ditto_colors %in% color_val)]
    }
    
    match_colors <- (coloring_algorithm != "none") # we want to match color if coloring_algorithm is not none
    
    if (!(coloring_algorithm %in% c("none", "left", "right", "advanced", graphing_columns))) {
        stop("Invalid coloring_algorithm. Options are 'none', 'left', 'right', 'advanced', or any value in graphing_columns.")
    }
    
    # warning if coloring_algorithm is not none but color_boxes is False
    if (!color_boxes && !color_bands && match_colors) {
        if (verbose) message("Warning: color_boxes and color_bands are False but coloring_algorithm is specified. boxes will not be colored.")
    }
    
    # Extract colors for each factor, assuming ditto_colors is long enough
    if (match_colors) {
        unused_colors <- ditto_colors
        first <- TRUE
        final_colors <- c()
        
        if (coloring_algorithm == "left") {
            n <- 1
            for (col_group in graphing_columns) {
                num_levels <- length(levels(df[[col_group]]))
                if (first) {
                    old_colors <- unused_colors[1:num_levels]
                    names(old_colors) <- 1:num_levels
                    unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                    final_colors <- c(final_colors, rev(old_colors))
                    first <- FALSE
                } else {
                    temp_colors <- find_group2_colors(df, ditto_colors, unused_colors, old_colors,
                                                      group1_name = paste0("col", n - 1, "_int"), group2_name = paste0("col", n, "_int"),
                                                      cutoff = cutoff
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                    old_colors <- temp_colors
                    names(old_colors) <- 1:num_levels
                    final_colors <- c(final_colors, rev(old_colors))
                }
                n <- n + 1
            }
        } else if (coloring_algorithm == "right") {
            n <- length(graphing_columns)
            for (col_group in rev(graphing_columns)) {
                num_levels <- length(levels(df[[col_group]]))
                if (first) {
                    old_colors <- unused_colors[1:num_levels]
                    names(old_colors) <- 1:num_levels
                    unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                    final_colors <- c(final_colors, rev(old_colors))
                    first <- FALSE
                } else {
                    temp_colors <- find_group2_colors(df, ditto_colors, unused_colors, old_colors,
                                                      group1_name = paste0("col", n, "_int"), group2_name = paste0("col", n - 1, "_int"),
                                                      cutoff = cutoff
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
                    old_colors <- temp_colors
                    names(old_colors) <- 1:num_levels
                    final_colors <- c(rev(old_colors), final_colors)
                    n <- n - 1
                }
            }
        } else if (coloring_algorithm == "advanced") {
            # check_python_setup_with_necessary_packages(necessary_packages_for_this_step = c("igraph", "leidenalg"), additional_message = "do not set coloring_algorithm to 'advanced'")
            final_colors <- find_colors_advanced(df, ditto_colors, graphing_columns, coloring_algorithm_advanced_option = coloring_algorithm_advanced_option, resolution = resolution)
        } else {
            col_group <- coloring_algorithm
            num_levels <- length(levels(df[[col_group]]))
            old_colors <- unused_colors[1:num_levels]
            names(old_colors) <- 1:num_levels
            unused_colors <- unused_colors[!(unused_colors %in% old_colors)]
            final_colors <- c(final_colors, rev(old_colors))
            ref_colors <- rev(old_colors)
            really_final_colors <- c()
            ref_group_n <- which(coloring_algorithm == graphing_columns)[[1]]
            
            for (col_group in graphing_columns) {
                if (!(col_group == coloring_algorithm)) {
                    col_group_n <- which(col_group == graphing_columns)[[1]]
                    num_levels <- length(levels(df[[col_group]]))
                    temp_colors <- find_group2_colors(df, ditto_colors, unused_colors, old_colors,
                                                      group1_name = paste0("col", ref_group_n, "_int"), group2_name = paste0("col", col_group_n, "_int"),
                                                      cutoff = cutoff
                    )
                    unused_colors <- unused_colors[!(unused_colors %in% temp_colors)]
                    final_colors <- c(final_colors, rev(temp_colors))
                    really_final_colors <- c(really_final_colors, rev(temp_colors))
                } else {
                    really_final_colors <- c(really_final_colors, ref_colors)
                }
            }
            final_colors <- really_final_colors
        }
    } else {
        remaining_colors <- ditto_colors
        final_colors <- c()
        for (col_group in graphing_columns) {
            num_levels <- length(levels(df[[col_group]]))
            if (length(remaining_colors) < num_levels) {
                if (verbose) message("Warning: Some colors will be recycled.")
                remaining_colors <- ditto_colors
            }
            old_colors <- remaining_colors[1:num_levels]
            final_colors <- c(final_colors, rev(old_colors))
            remaining_colors <- remaining_colors[(1 + num_levels):length(remaining_colors)]
        }
    }
    
    if (!(is.null(color_val))) {
        final_value_order <- c()
        for (col_int in seq_along(graphing_columns)) {
            int_name <- paste0("col", col_int, "_int")
            group_name <- graphing_columns[[col_int]]
            
            curr_label <- as.character(unique(df[order(df[[int_name]]), ][[group_name]]))
            
            final_value_order <- c(final_value_order, rev(curr_label))
        }
        names(final_colors) <- final_value_order
        
        for (box_val in names(color_val)) {
            if (box_val %in% names(final_colors)) {
                box_val_color <- color_val[names(color_val) == box_val][1]
                # find where value is in final colors
                val_match <- which(box_val == names(final_colors))
                
                # JR added
                to_change <- c()
                for (vm in val_match) {
                    matched_color <- final_colors[vm]
                    maybe_to_change <- which(matched_color == final_colors)
                    
                    for (matched in maybe_to_change) {
                        testing <- names(final_colors)[matched]
                        if (!(testing %in% names(color_val)) || (testing == box_val)) {
                            to_change <- c(to_change, matched)
                        }
                    }
                }
                # JR added
                
                # JR commented out
                # matched_color <- final_colors[val_match[1]]
                # maybe_to_change <- which(matched_color == final_colors)
                # to_change <- c()
                # for (matched in maybe_to_change) {
                #     testing <- names(final_colors)[matched]
                #     if (!(testing %in% names(color_val)) | (testing == box_val)) {
                #         to_change <- c(to_change, matched)
                #     }
                # }
                # JR commented out
                
                final_colors[as.integer(to_change)] <- box_val_color
            }
        }
    }
    
    remaining_colors <- ditto_colors[!(ditto_colors %in% final_colors)]
    
    final_colors_legend <- final_colors
    # generate label names
    final_label_names <- c()
    for (col_int in seq_along(graphing_columns)) {
        int_name <- paste0("col", col_int, "_int")
        group_name <- graphing_columns[[col_int]]
        
        curr_label <- as.character(unique(df[order(df[[int_name]]), ][[group_name]]))
        
        final_label_names <- c(final_label_names, rev(curr_label))
    }
    names(final_colors_legend) <- final_label_names
    
    # if color_bands, add to named list
    
    if (!is.null(color_band_list)) {
        final_colors_legend <- c(color_band_list, final_colors_legend)
    } else {
        if (!is.null(color_band_column)) {
            if (!(color_band_column %in% graphing_columns)) {
                num_levels <- length(unique(df[[color_band_column]]))
                color_band_list <- remaining_colors[1:num_levels]
                names(color_band_list) <- unique(df[[color_band_column]])
                final_colors_legend <- c(color_band_list, final_colors_legend)
            }
        }
    }
    # remove duplicate names
    final_colors_legend <- final_colors_legend[!duplicated(names(final_colors_legend))]
    # remove duplicate dims
    temp_df <- df # [1:as.integer(dim(df)[1]/2),1:dim(df)[2]]
    
    # uncomment to attempt mapping
    p <- ggplot(data = temp_df, aes(y = value), )
    for (x in seq_along(graphing_columns)) {
        p$mapping[[paste0("axis", x)]] <- sym(paste0("col", x, "_int"))
    }
    
    if (color_bands) {
        if (is.null(color_band_column)) {
            color_band_column <- graphing_columns[1]
        }
        if (color_band_boundary) {
            p <- p +
                geom_alluvium(aes(fill = !!sym(color_band_column), color = !!sym(color_band_column)),
                              alpha = alluvial_alpha
                )
        } else {
            p <- p +
                geom_alluvium(aes(fill = !!sym(color_band_column)), alpha = alluvial_alpha)
        }
        p <- p + scale_fill_manual(values = final_colors_legend)
    } else {
        if (color_band_boundary) {
            p <- p + geom_alluvium(color = "grey2", alpha = alluvial_alpha)
        } else {
            p <- p + geom_alluvium(alpha = alluvial_alpha)
        }
    }
    
    if (color_band_boundary) {
        if (add_legend) {
            p <- p + scale_color_manual(values = final_colors_legend, guide = "none")
        } else {
            p <- p + scale_color_manual(values = color_band_list, guide = "none")
        }
    }
    
    
    if (color_boxes) {
        if (add_legend) {
            p <- p + geom_stratum(width = box_width, aes(fill = after_stat(!!sym("final_label_names"))), linewidth = box_line_width) + scale_fill_manual(values = final_colors_legend)
        } else {
            p <- p + geom_stratum(width = box_width, fill = final_colors, linewidth = box_line_width)
        }
    } else {
        p <- p + geom_stratum(width = box_width, linewidth = box_line_width)
    }
    
    
    if (!(include_labels_in_boxes == FALSE)) {
        if (auto_adjust_text) {
            p <- p +
                geom_fit_text(
                    reflow = TRUE, stat = StatStratum, width = text_width, min.size = min_text, size = text_size,
                    aes(label = after_stat(final_label_names))
                )
        } else {
            p <- p +
                geom_text(stat = StatStratum, aes(label = after_stat(final_label_names)), size = text_size)
        }
    }
    
    top_y <- 0
    for (test_x in unique(df$x)) {
        curr_y <- df %>%
            filter(x == test_x) %>%
            group_by(y) %>%
            summarise(total = sum(value), .groups = "drop") %>%
            arrange(desc(total)) %>%
            mutate(cum_y = cumsum(total)) %>%
            pull(cum_y) %>%
            max()
        top_y <- max(curr_y, top_y)
    } # top_y1 and top_y2 are probably the same
    
    if (include_axis_titles) {
        p <- p +
            scale_x_continuous(
                breaks = 1:length(graphing_columns), labels = graphing_columns,
                position = "top"
            )
    }
    
    if (include_group_sizes) {
        offset_below <- top_y * 0.075
        x <- 1
        for (col_group in graphing_columns) {
            p <- p +
                annotate("text", x = x, y = -offset_below, label = length(levels(df[[col_group]])), hjust = 0.5, size = 5) # Adjust x, y for Scanpy
            x <- x + 1
        }
    }
    
    if (flip_xy) {
        p <- p + ggplot2::coord_flip()
    }
    
    p <- p +
        theme_void() +
        theme(
            text = element_text(family = "sans"))
    
    if (add_legend) {
        p <- p + theme(legend.position = legend_loc) + labs(fill = "")
    } else {
        p <- p + theme(legend.position = "none")
    }
    
    if (keep_x_labels) {
        p <- p + theme(axis.text.x = element_text(size = axis_text_size, vjust = axis_text_vjust)) # vjust adjusts the vertical position of column titles)
    } 
    if (keep_y_labels) {
        p <- p + theme(axis.text.y = element_text(size = axis_text_size, vjust = axis_text_vjust))
    }
    
    
    if (!is.null(output_plot_path)) {
        if (verbose) message(sprintf("Saving plot to: %s", output_plot_path))
        ggsave(output_plot_path,
               plot = p,
               height = save_height, width = save_width, dpi = dpi, bg = "white"
        )
    }
    
    return(p)
}
