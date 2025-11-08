test_that("plot_alluvial returns a ggplot object", {
    set.seed(42)
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    p <- plot_alluvial(df, sorting_algorithm="tsp", coloring_algorithm="left")
    expect_s3_class(p, "ggplot")
})

test_that("plot_alluvial returns an error with 0 rows", {
    df <- data.frame()
    expect_error(plot_alluvial(df, sorting_algorithm="tsp"), "no rows")
})

test_that("plot_alluvial returns an error with 1 column", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, sorting_algorithm="tsp"), "at least 2 columns")
})

test_that("plot_alluvial returns an error with 2 columns, where column1 is not in df", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, graphing_columns = c("bad_col", "method2"), sorting_algorithm="tsp"), "column 'bad_col' is not a column in the dataframe")
})

test_that("plot_alluvial returns an error with 2 columns, where column2 is not in df", {
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_error(plot_alluvial(df, graphing_columns = c("method1", "bad_col"), sorting_algorithm="tsp"), "column 'bad_col' is not a column in the dataframe")
})

test_that("plot_alluvial returns an error with 3 columns, where graphing_columns is NULL", {
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE),
        method3 = sample(1:3, 100, TRUE)
    )
    expect_error(plot_alluvial(df, graphing_columns = NULL, sorting_algorithm="tsp"), "dataframe has more than 2 columns")
})

test_that("plot_alluvial returns an error with 3 columns, where column2 is not in df", {
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE),
        method3 = sample(1:3, 100, TRUE)
    )
    expect_error(plot_alluvial(df, graphing_columns = c("method1", "bad_col"), sorting_algorithm="tsp"), "column 'bad_col' is not a column in the dataframe")
})

test_that("plot_alluvial works with a df with 2 columns, both column1 and column2 NULL", {
    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    expect_s3_class(plot_alluvial(df, sorting_algorithm="tsp", coloring_algorithm="left"), "ggplot")
})

test_that("plot_alluvial works with 4 columns, column1 and column2 specified", {
    set.seed(42)
    df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE),
        meta = rnorm(100),
        sample_id = paste0("sample_", 1:100)
    )
    expect_s3_class(plot_alluvial(df, graphing_columns = c("method1", "method2"), sorting_algorithm="tsp", coloring_algorithm="left"), "ggplot")
})

test_that("plot_alluvial works with a pre-aggregated df with weight column", {
    set.seed(42)
    # Generate raw data
    raw_df <- data.frame(
        method1 = sample(1:3, 100, TRUE),
        method2 = sample(1:3, 100, TRUE)
    )

    # Aggregate by combination
    df <- as.data.frame(dplyr::count(raw_df, method1, method2, name = "value"))

    expect_s3_class(plot_alluvial(df, graphing_columns = c("method1", "method2"), column_weights = "value", sorting_algorithm="tsp", coloring_algorithm="left"), "ggplot")
})

test_that("VDIFFR - plot_alluvial works with 2 columns, unsorted", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    p <- plot_alluvial(df, graphing_columns = c("method1", "method2"), sorting_algorithm = "none", coloring_algorithm="left")
    vdiffr::expect_doppelganger("basic_alluvial_plot_unsorted", p)
})

test_that("VDIFFR - plot_alluvial works with 2 columns, greedy_wolf", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    p <- plot_alluvial(df, graphing_columns = c("method1", "method2"), sorting_algorithm = "greedy_wolf", coloring_algorithm="left")
    vdiffr::expect_doppelganger("basic_alluvial_plot_WOLF", p)
})

test_that("VDIFFR - plot_alluvial works with 2 columns, greedy_wblf", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    p <- plot_alluvial(df, graphing_columns = c("method1", "method2"), sorting_algorithm = "greedy_wblf", coloring_algorithm="left")
    vdiffr::expect_doppelganger("basic_alluvial_plot_WBLF", p)
})

test_that("VDIFFR - plot_alluvial works with 2 columns, TSP", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")
    
    set.seed(42)
    df <- data.frame(method1 = sample(1:3, 100, TRUE), method2 = sample(1:3, 100, TRUE))
    p <- plot_alluvial(df, graphing_columns = c("method1", "method2"), sorting_algorithm = "tsp")
    vdiffr::expect_doppelganger("basic_alluvial_plot_TSP", p)
})


make_more_tsp_2_layer_df <- function() {
    df <- data.frame(
        tissue = c(
            1, 1, 1,
            2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3,
            4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5
        ),
        cluster = c(
            6, 6, 7,
            6, 7, 7, 7, 7, 7,
            6, 8, 8, 8, 8, 8, 8,
            8, 8,
            8, 8, 8, 8, 8, 8, 8, 8, 8
        )
    )
    column1 <- "tissue"
    column2 <- "cluster"
    graphing_columns <- c(column1, column2)

    list(
        df = df,
        graphing_columns = graphing_columns
    )
}

test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 2 layers, unsorted", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_2_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "none", coloring_algorithm="left", color_bands = TRUE)

    vdiffr::expect_doppelganger("more_tsp_2layer_unsorted", p)
})


test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 2 layers, tsp, optimize_column_order_per_cycle FALSE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_2_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "tsp", color_bands = TRUE, optimize_column_order_per_cycle = FALSE)

    vdiffr::expect_doppelganger("more_tsp_2layer_tsp_oFALSE", p)
})

test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 2 layers, tsp, optimize_column_order_per_cycle TRUE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_2_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "tsp", color_bands = TRUE, optimize_column_order_per_cycle = TRUE)

    vdiffr::expect_doppelganger("more_tsp_2layer_tsp_oTRUE", p)
})




make_more_tsp_3_layer_df <- function() {
    df <- data.frame(
        tissue = c(
            "BRAIN", "BRAIN", "BRAIN",
            "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH",
            "HEART", "HEART", "HEART", "HEART", "HEART", "HEART", "HEART",
            "T CELL", "T CELL",
            "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL"
        ),
        cluster = c(
            1, 1, 2,
            1, 2, 2, 2, 2, 2,
            1, 3, 3, 3, 3, 3, 3,
            4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4
        ),
        sex = c(
            "male", "female", "male",
            "female", "male", "female", "female", "male", "female",
            "male", "female", "male", "female", "male", "female", "male",
            "female", "male",
            "male", "male", "male", "male", "male", "male", "male", "male", "male"
        )
    )
    graphing_columns <- c("tissue", "cluster", "sex")

    list(
        df = df,
        graphing_columns = graphing_columns
    )
}

make_more_tsp_3_layer_df_with_2_identical_layers <- function() {
    df <- data.frame(
        tissue = c(
            "BRAIN", "BRAIN", "BRAIN",
            "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH",
            "HEART", "HEART", "HEART", "HEART", "HEART", "HEART", "HEART",
            "T CELL", "T CELL",
            "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL"
        ),
        sex = c(
            "male", "female", "male",
            "female", "male", "female", "female", "male", "female",
            "male", "female", "male", "female", "male", "female", "male",
            "female", "male",
            "female", "male", "female", "male", "female", "male", "female", "female", "male"
        ),
        cluster = c(
            "BRAIN", "BRAIN", "BRAIN",
            "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH", "STOMACH",
            "HEART", "HEART", "HEART", "HEART", "HEART", "HEART", "HEART",
            "T CELL", "T CELL",
            "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL", "B CELL"
        )
    )
    graphing_columns <- c("tissue", "cluster", "sex")

    list(
        df = df,
        graphing_columns = graphing_columns
    )
}

test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 3 layers, unsorted", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "none", coloring_algorithm="left", color_bands = TRUE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_tsp_3layer_unsorted", p)
})



test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 3 layers, tsp, optimize_column_order FALSE, optimize_column_order_per_cycle FALSE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "tsp", color_bands = TRUE, optimize_column_order = FALSE, optimize_column_order_per_cycle = FALSE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_tsp_3layer_unsorted_oFALSE_ocFALSE", p)
})


test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 3 layers, tsp, optimize_column_order FALSE, optimize_column_order_per_cycle TRUE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "tsp", color_bands = TRUE, optimize_column_order = FALSE, optimize_column_order_per_cycle = TRUE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_tsp_3layer_unsorted_oFALSE_ocTRUE", p)
})

test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 3 layers, tsp, optimize_column_order TRUE, optimize_column_order_per_cycle FALSE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "tsp", color_bands = TRUE, optimize_column_order = TRUE, optimize_column_order_per_cycle = FALSE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_tsp_3layer_unsorted_oTRUE_ocFALSE", p)
})

test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 3 layers, tsp, optimize_column_order TRUE, optimize_column_order_per_cycle TRUE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_3_layer_df()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "tsp", color_bands = TRUE, optimize_column_order = TRUE, optimize_column_order_per_cycle = TRUE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_tsp_3layer_unsorted_ooTRUE_ocTRUE", p)
})


test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 3 layers with 2 identical layers, unsorted", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "none", color_bands = TRUE)

    vdiffr::expect_doppelganger("more_tsp_3layer_2ident", p)
})

test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 3 layers with 2 identical layers, tsp, optimize_column_order FALSE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "tsp", color_bands = TRUE, optimize_column_order = FALSE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_tsp_3layer_with_2ident_oFALSE", p)
})


test_that("VDIFFR - plot_alluvial, more_tsp.Rmd, 3 layers with 2 identical layers, tsp, optimize_column_order TRUE", {
    skip_if_not(requireNamespace("vdiffr", quietly = TRUE), "vdiffr not installed")

    set.seed(42)
    input <- make_more_tsp_3_layer_df_with_2_identical_layers()
    df <- input$df
    graphing_columns <- input$graphing_columns

    p <- plot_alluvial(df, graphing_columns = graphing_columns, sorting_algorithm = "tsp", color_bands = TRUE, optimize_column_order = TRUE, weight_scalar = 1)

    vdiffr::expect_doppelganger("more_tsp_3layer_with_2_ident_oTRUE", p)
})
