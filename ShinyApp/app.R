library(shiny)
library(shinydashboard)
library(magrittr)
library(ggplot2)
library(shinyWidgets)
library(shinyBS)
library(dplyr)
library(Seurat)
library(ggtext)

source(here::here("plot_theme.R"))

# once a umap is plotted, it will be cached for speedup.
# We have to limit the max RAM this is using.
shinyOptions(cache = cachem::cache_mem(max_size = 200e6))


# load the srt object.
# reactivity not really important here, because it is used all the time.
# this is the normal seurat object, but all genes that are not expressed were removed and also
# all dimensionality reductions to reduce the size and speedup loading
# see "filtering_steps.R" for details.
seurat <- readRDS(here::here("data", "filtered_seurat_object.rds"))
seurat_RNAi <- readRDS(here::here("data", "filtered_seurat_object_RNAi.rds"))

# see "create_feature_mapping" on details how this was created.
feature_mapping <- read.csv(here::here("data", "fbid_gene_mapping.csv"))


#############
# FUNCTIONS #
#############

# Helper function to map gene names
map_gene_name <- function(gene_name, mapping) {
    if (grepl("FBgn\\d+", gene_name)) {
        return(mapping %>%
            filter(fbid == gene_name) %>%
            pull(symbol))
    }
    return(gene_name)
}

# Helper function to get condition-specific title and styling
get_condition_config <- function(split_conditions, perturbation_count) {
    config <- list(
        title = "Combined Dataset",
        scale_transform = NULL,
        title_style = NULL
    )

    switch(split_conditions,
        "combined" = {
            config$title <- "Combined Dataset"
            if (perturbation_count == 2) config$scale_transform <- "scale_x_reverse"
        },
        "ctrl" = {
            config$title <- "Control Condition"
            config$title_style <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
            if (perturbation_count == 2) config$scale_transform <- "scale_x_reverse"
        },
        "notch" = {
            config$title <- "*Notch*<sup><i>sgRNAx2</i></sup>"
            config$title_style <- theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
            config$scale_transform <- "scale_x_reverse"
        },
        "NotchRNAi" = {
            config$title <- "*Notch*<sup><i>RNAi</i></sup>"
            config$title_style <- theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
        },
        "NotchCphRNAi" = {
            config$title <- "*Cph*<sup><i>RNAi</i></sup>+*Notch*<sup><i>RNAi</i></sup>"
            config$title_style <- theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
        },
        "CphUp" = {
            config$title <- "*Cph*<sup>up</sup>"
            config$title_style <- theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
        }
    )

    return(config)
}

# Helper function to apply scale transformations
apply_scale_transform <- function(plot, transform) {
    if (is.null(transform)) return(plot)

    switch(transform,
        "scale_x_reverse" = plot + scale_x_reverse(),
        "scale_y_reverse" = plot + scale_y_reverse(),
        plot
    )
}

# this function returns the actual ggplot
annotation_plot_call <- function(data, dimensions) {
    ggplot(data, aes(x = !!sym(dimensions[1]), y = !!sym(dimensions[2]), color = celltype_manual)) +
        geom_point(shape = ".") +
        small_axis(substr(dimensions[1], 1, nchar(dimensions[1]) - 2), fontsize = 8, arrow_length = 20) +
        theme(
            legend.text = element_text(size = 13.5),
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.4, "cm"),
            legend.title = element_blank()
        ) +
        scale_color_manual(values = celltype_colors) +
        guides(color = guide_legend(override.aes = list(size = 3, shape = 19)))
}

# this function handles different cases and then calls `annotation_plot_call`.
annotation_plot_fun <- function(srt, dim_red, split_conditions) {
    if (dim_red == "PCA") {
        dim_red <- "PC"
    }

    dimensions <- c(paste0(dim_red, "_1"), paste0(dim_red, "_2"))
    perturbation_count <- n_distinct(srt$perturbation)

    # generate the plot data
    plot_data <- srt@meta.data[, c("celltype_manual", dimensions, "perturbation")] %>%
        dplyr::sample_frac()

    # Filter data based on condition
    if (split_conditions != "combined") {
        plot_data <- plot_data %>% filter(perturbation == split_conditions)
    }

    # Get condition configuration
    config <- get_condition_config(split_conditions, perturbation_count)

    # Create base plot
    p <- annotation_plot_call(plot_data, dimensions) +
        ggtitle(config$title)

    # Apply title styling if specified
    if (!is.null(config$title_style)) {
        p <- p + config$title_style
    }

    # Apply scale transformation if specified
    p <- apply_scale_transform(p, config$scale_transform)

    return(p)
}

# this function returns the actual ggplot
feature_plot_call <- function(data, dimensions, gene) {
    ggplot(data, aes(x = !!sym(dimensions[1]), y = !!sym(dimensions[2]), color = !!sym(gene))) +
        geom_point(shape = ".") +
        small_axis(substr(dimensions[1], 1, nchar(dimensions[1]) - 2), fontsize = 8, arrow_length = 20) +
        theme(
            legend.key.width = unit(0.9, "line"),
            legend.text = element_text(size = 9),
            legend.key.height = unit(0.7, "line")
        ) +
        scale_color_gradient(
            low = "lightgrey", high = "darkred",
            name = "Expression\n(logcounts)\n"
        )
}

# this function handles different cases and then calls `feature_plot_call`.
feature_plot_fun <- function(srt, gene, dim_red, order_cells, split_conditions) {
    if (dim_red == "PCA") {
        dim_red <- "PC"
    }
    dimensions <- c(paste0(dim_red, "_1"), paste0(dim_red, "_2"))
    perturbation_count <- n_distinct(srt$perturbation)

    # Prepare plot data
    mdata <- srt@meta.data[, c("celltype_manual", dimensions, "perturbation")]
    exp_data <- FetchData(srt, gene)
    colnames(exp_data) <- gene
    stopifnot(all(rownames(mdata) == rownames(exp_data)))

    plot_data <- cbind(mdata, exp_data)

    # Order or shuffle data
    if (order_cells) {
        plot_data <- plot_data %>%
            dplyr::arrange(., !!sym(gene))
    } else {
        plot_data <- plot_data %>%
            dplyr::sample_frac()
    }

    # Filter data based on condition
    if (split_conditions != "combined") {
        # Special case: CphUp condition uses NotchCphRNAi data
        filter_condition <- if (split_conditions == "CphUp") "NotchCphRNAi" else split_conditions
        plot_data <- plot_data %>% filter(perturbation == filter_condition)
    }

    # Get condition configuration
    config <- get_condition_config(split_conditions, perturbation_count)

    # Create base plot
    p <- feature_plot_call(plot_data, dimensions, gene) +
        ggtitle(config$title)

    # Apply title styling if specified
    if (!is.null(config$title_style)) {
        p <- p + config$title_style
    }

    # Apply scale transformation if specified
    p <- apply_scale_transform(p, config$scale_transform)

    return(p)
}

jitter_plot_fun <- function(srt, gene, celltypes, summary_stats, split_conditions, conditions = NULL, use_rnai_colors = FALSE) {
  mdata <- srt@meta.data[, c("celltype_manual", "perturbation")]
  exp_data <- FetchData(srt, gene)
  colnames(exp_data) <- gene
  stopifnot(all(rownames(mdata) == rownames(exp_data)))

  plot_data <- cbind(mdata, exp_data) %>%
    filter(celltype_manual %in% celltypes)

  # Filter by specific conditions if provided (for RNAi plots)
  if (!is.null(conditions)) {
    plot_data <- plot_data %>% filter(perturbation %in% conditions)
  }

  # Helper function to add summary statistics
  add_summary_stats <- function(plot, stats, grouped = FALSE) {
    if (is.null(stats)) return(plot)

    if ("Include mean" %in% stats) {
      if (grouped) {
        plot <- plot + stat_summary(
          fun = "mean", geom = "point", aes(group = perturbation),
          color = "black", size = 3, shape = 8, stroke = 1.2,
          position = position_dodge2(width = 0.75)
        )
      } else {
        plot <- plot + stat_summary(
          fun = "mean", geom = "point",
          color = "black", size = 3, shape = 8, stroke = 1.2
        )
      }
    }

    if ("Include median" %in% stats) {
      if (grouped) {
        plot <- plot + stat_summary(
          fun = "median", geom = "point", aes(group = perturbation),
          color = "black", size = 3, shape = 2, stroke = 1.2,
          position = position_dodge2(width = 0.75)
        )
      } else {
        plot <- plot + stat_summary(
          fun = "median", geom = "point",
          color = "black", size = 3, shape = 2, stroke = 1.2
        )
      }
    }

    return(plot)
  }

  if (split_conditions == "combined") {
    p <- ggplot(plot_data, aes(x = celltype_manual, y = !!sym(gene), color = celltype_manual, fill = celltype_manual)) +
      geom_jitter(size = 0.65, stroke = 0.6, width = 0.2, shape = 21) +
      geom_violin(
        adjust = 1, trim = TRUE, color = "black", show.legend = F,
        scale = "width", fill = "#00000000",
        width = 0.72
      ) +
      theme_Publication() +
      theme(
        axis.title.x = element_blank(),
        legend.position = "none"
      ) +
      ylab("Expression\n(logcounts)") +
      scale_color_manual(values = celltype_colors) +
      scale_fill_manual(values = alpha(celltype_colors, 0.4)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08)))

    p <- add_summary_stats(p, summary_stats, grouped = FALSE)

  } else if (split_conditions == "split" || !is.null(conditions)) {
    p <- ggplot(plot_data, aes(x = celltype_manual, y = !!sym(gene))) +
      geom_point(
        mapping = aes(color = perturbation, fill = perturbation), size = 0.65, stroke = 0.4, shape = 21,
        position = position_jitterdodge()
      ) +
      geom_violin(
        mapping = aes(fill = perturbation),
        adjust = 1, trim = TRUE,
        scale = "width", show.legend = F,
        width = 0.75
      ) +
      theme_Publication_side_legend() %+%
      theme(
        axis.title.x = element_blank(),
        legend.text = element_markdown(size = 15),
        legend.title = element_text(size = 16)
      ) +
      ylab("Expression\n(logcounts)") +
      guides(
        color = guide_legend(
          title = "Condition",
          override.aes = list(size = 2, alpha = 1, shape = NULL),
          hjust = 0
        ),
        fill = "none"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08)))

    p <- add_summary_stats(p, summary_stats, grouped = TRUE)

    # Apply color schemes
    if (use_rnai_colors) {
      p <- p +
        scale_color_manual(values = color_mapping_rnai) +
        scale_fill_manual(values = alpha(color_mapping_rnai, 0))
    } else if (n_distinct(srt$perturbation) == 3) {
      p <- p +
        scale_color_manual(
          values = c("ctrl" = color_mapping[1], "NotchRNAi" = color_mapping[4], NotchCphRNAi = "#00ba38"),
          labels = c(
            ctrl = "Control",
            NotchRNAi = "*Notch*<sup><i>RNAi</i></sup>",
            NotchCphRNAi = "*Cph*<sup><i>RNAi</i></sup>+*Notch*<sup><i>RNAi</i></sup>"
          )
        ) +
        scale_fill_manual(values = alpha(c(color_mapping[1], color_mapping[4], "#00ba38"), 0))
    } else if (n_distinct(srt$perturbation) == 2) {
      p <- p +
        scale_color_manual(
          values = c("ctrl" = color_mapping[1], "notch" = color_mapping[4]),
          labels = c(
            ctrl = "Control",
            notch = "*Notch*<sup><i>sgRNAx2</i></sup>"
          )
        ) +
        scale_fill_manual(values = alpha(c("#00000000", "#00000000"), 0))
    }
  }

  return(p + ggtitle(gene))
}

# Wrapper function for RNAi-specific jitter plots (for backward compatibility)
jitter_plot_fun_rnai <- function(srt, gene, celltypes, summary_stats, conditions) {
  return(jitter_plot_fun(srt, gene, celltypes, summary_stats, "split", conditions, use_rnai_colors = TRUE))
}

vln_gene_by_region <- function(obj, gene, slot = "data") {
  df <- FetchData(obj, vars = c(gene, "region_prediction", "perturbation"), slot = slot)
  names(df)[1] <- "expr"

  df <- df[df$perturbation == "ctrl" &
             !is.na(df$region_prediction) &
             df$region_prediction != "region prediction not available", ]

  levs <- unique(df$region_prediction)
  if (all(grepl("^R\\d+$", levs))) {
    levs <- levs[order(as.integer(sub("^R", "", levs)))]
  } else {
    levs <- sort(levs)
  }
  df$region_prediction <- factor(df$region_prediction, levels = levs)

  strong_colors <- c(
    "R1" = "#D55E00",
    "R2" = "#E69F00",
    "R3" = "#009E73",
    "R4" = "#0072B2",
    "R5" = "#CC79A7"
  )

  ggplot(df, aes(x = region_prediction, y = expr, color = region_prediction, fill = region_prediction)) +
    geom_jitter(size = 0.65, stroke = 0.6, width = 0.2, shape = 21) +
    geom_violin(
      adjust = 1, trim = TRUE, color = "black", show.legend = F,
      scale = "width", fill = "#00000000",
      width = 0.72
    ) +
    theme_Publication() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ylab("Expression\n(logcounts)") +
    scale_color_manual(values = strong_colors) +
    scale_fill_manual(values = alpha(strong_colors, 0.4)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    ggtitle(gene)
}

ui <- dashboardPage(
    dashboardHeader(
        title = "", titleWidth = 115,
        # Set height of dashboardHeader
        tags$li(
            class = "dropdown",
            tags$style(".main-header {max-height: 20px}"),
            tags$style(".main-header .logo {height: 20px;}"),
            tags$style(".sidebar-toggle {height: 20px; padding-top: 1px !important;}"),
            tags$style(".navbar {min-height:20px !important}")
        )
    ),
    dashboardSidebar(
        tags$style(".left-side, .main-sidebar {padding-top: 20px}"),
        sidebarMenu(
            menuItem("About", tabName = "about", icon = icon("fas fa-circle-info", class = "fa-lg")),
            menuItem(HTML("<i>Notch<sup>sgRNAx2</i></sup><br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dataset"),
                tabName = "dataset_notchKO", icon = icon("fas fa-magnifying-glass", class = "fa-lg")
            ),
            menuItem(HTML("RNAi<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;dataset"),
                tabName = "dataset_RNAi", icon = icon("fas fa-magnifying-glass", class = "fa-lg")
            ),
            menuItem(HTML("Regional<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Expression"),
                tabName = "regional_expression", icon = icon("fas fa-chart-area", class = "fa-lg")
            ),
            menuItem("Contact", tabName = "contact", icon = icon("far fa-address-book", class = "fa-lg"))
        ),
        width = 115
    ),
    dashboardBody(
        tags$style(
            HTML("
           .col-sm-4 {
             padding: 5px !important; /* Adjust padding value as needed */
           }
           .col-sm-3 {
             padding: 5px !important; /* Adjust padding value as needed */
           }
           .col-sm-8 {
             padding: 5px !important; /* Adjust padding value as needed */
           }
           .box-title {

             font-size: 14px !important;
             font-weight: bold !important; /* Add this line to make title bold */

           }
           .box-header {
             padding: 3px 8px !important; /* Adjust padding values as needed */
           }
           .box.box-solid.box-primary {
             margin-bottom: 0px !important; /* Adjust margin value as needed */
           }
           .irs-handle.single {
           width: 5px;
           height: 5px;
           }
           /* Reduce margin and padding of input elements inside the box */
          .box .form-group {
            margin-top: 1px !important;
            margin-bottom: 1px !important;
            padding-top: 1px !important;
            padding-bottom: 1px !important;
          }
          /* Further adjustments for radio buttons */
          .box .radio {
            margin-top: 1px !important;
            margin-bottom: 1px !important;
          }
          /* Adjust checkbox group input */
          .box .checkbox {
            margin-top: 1px !important;
            margin-bottom: 1px !important;
          }
          /* Action button */
          .box .btn {
            margin-top: 10px !important;  /* Give the button some top margin */
            margin-bottom: 1px !important;
          }

          .box .selectize-control {
            margin-top: 0px !important;
            margin-bottom: 0px !important;
            padding-top: 0px !important;
            padding-bottom: 0px !important;
          }
           ")
        ),
        tabItems(
            tabItem(
                tabName = "about",

            h4(
                "An intrinsic autoinhibitory feedback loop preserves intestinal stem cell maintenance and fate commitment",
                align = "center",
                style = "width: 80%; margin: 0 auto; font-weight: bold; font-style: italic;"
            ),

            div(style = "margin-bottom: 10px;"),

            # Authors (numbers = affiliations; *, †, # as in the Word doc)
            h5(
                HTML(paste0(
                "Siamak Redhai<sup>1,2,3</sup><sup>*</sup><sup>#</sup>, ",
                "Nick Hirschmüller<sup>5</sup><sup>*</sup><sup>†</sup>, ",
                "Tianyu Wang<sup>1,2,3</sup>, ",
                "Tyler Jackson<sup>6</sup>, ",
                "Stefan Peidli<sup>5</sup>, ",
                "Erica Valentani<sup>1,2,3</sup>, ",
                "Shivohum Bahuguna<sup>1,2,3</sup>, ",
                "Svenja Leible<sup>1,2,3</sup>, ",
                "Sarina Möller<sup>1,3</sup>, ",
                "Christina Ko<sup>8</sup>, ",
                "Michaela Holzem<sup>1,2,3</sup>, ",
                "Sviatoslav Kharuk<sup>5,6</sup>, ",
                "Lea Bräckow<sup>1,3</sup>, ",
                "Fillip Port<sup>1,2,3</sup>, ",
                "David Ibberson<sup>4</sup>, ",
                "Hongji Le<sup>7</sup>, ",
                "Wolfgang Huber<sup>5</sup><sup>#</sup>, ",
                "Michael Boutros<sup>1,2,3</sup><sup>#</sup>"
                )),
                align = "center",
                style = "max-width: 1000px; margin: 0 auto; color: #555;"
            ),

            div(style = "height: 6px;"),

            # Footnotes (centered)
            h6(
                HTML(paste0(
                "<em>*</em> Contributed equally",
                "<br><strong>#</strong> Corresponding authors: ",
                "<a href='mailto:siamak.redhai@dkfz.de'>siamak.redhai@dkfz.de</a>, ",
                "<a href='mailto:wolfgang.huber@embl.org'>wolfgang.huber@embl.org</a>, ",
                "<a href='mailto:m.boutros@dkfz.de'>m.boutros@dkfz.de</a>"
                )),
                align = "center",
                style = "width: 90%; margin: 0 auto; color: #666;"
            ),

            div(style = "height: 8px;"),

            # Affiliations (2 columns, 4 affiliations each)
            div(
                HTML(paste0(
                "<div style='max-width: 920px; margin: 0 auto; color:#444; font-size: 0.9em; ",
                "display: flex; justify-content: space-between; gap: 20px;'>",
                "<div style='flex: 1;'>",
                "<p style='margin: 3px 0;'><sup>1</sup> German Cancer Research Center (DKFZ), Division Signaling and Functional Genomics, ",
                "Im Neuenheimer Feld 580, 69120 Heidelberg, Germany</p>",
                "<p style='margin: 3px 0;'><sup>2</sup> Heidelberg University, Institute of Human Genetics, Medical Faculty Heidelberg, ",
                "Im Neuenheimer Feld 366, 69120 Heidelberg, Germany</p>",
                "<p style='margin: 3px 0;'><sup>3</sup> Heidelberg University, Department of Cell and Molecular Biology, Medical Faculty Mannheim ",
                "&amp; BioQuant, Im Neuenheimer Feld 267, 69120 Heidelberg, Germany</p>",
                "<p style='margin: 3px 0;'><sup>4</sup> CellNetworks Core Technology Platform (CCTP), Deep Sequencing Labor, ",
                "Im Neuenheimer Feld 267, 69120 Heidelberg, Germany</p>",
                "</div>",
                "<div style='flex: 1;'>",
                "<p style='margin: 3px 0;'><sup>5</sup> European Molecular Biology Laboratory (EMBL), Meyerhofstraße 1, 69117 Heidelberg, Germany</p>",
                "<p style='margin: 3px 0;'><sup>6</sup> Department of Biochemistry and Biotechnology, Vasyl Stefanyk Precarpathian National University, ",
                "Ivano-Frankivsk, Ukraine</p>",
                "<p style='margin: 3px 0;'><sup>7</sup> Huffington Center on Aging &amp; Department of Molecular and Human Genetics, ",
                "Baylor College of Medicine, Houston, TX 77030, USA</p>",
                "<p style='margin: 3px 0;'><sup>8</sup> Rice University, BioSciences Department, Houston, TX 77030, USA</p>",
                "</div>",
                "</div>"
                )),
                align = "center"
            ),

            div(style = "margin-bottom: 30px;"),

            h3("Abstract", style = "font-weight: bold"),
            p("Intestinal stem cells (ISCs) continuously renew the gut epithelium by generating regionally specialized cell types, but how they commit to distinct lineages remains poorly defined. Here, we identify a self-limiting transcriptional program, mediated by the zinc-finger transcription factor Chronophage (Cph), that drives ISC maintenance and differentiation into enteroendocrine (EE) cells in different regions of the Drosophila midgut. We show that Cph expression is transiently induced by the proneural factor scute at the onset of ISC-to-EE fate specification. Genetic and single-cell transcriptomics approaches revealed that Cph is required to intrinsically remodel the transcriptome of ISCs and sustains normal lifespan. Genome-wide chromatin profiling demonstrated that Cph directly binds and regulates the expression of key proliferation and differentiation genes, while simultaneously repressing its own expression. This autoinhibitory feedback safeguards ISCs from undergoing autophagy and cell death, thus ensuring proliferation and differentiation are faithfully executed. Our findings highlight a key mechanism that balances ISC maintenance and lineage commitment."),
            p(
                "This Shiny App allows users to interactively explore the datasets accompanying the publication. ",
                "The source code for the app can be found on ",
                a("GitHub", href = "https://github.com/nickhir/Chronophage/tree/main/ShinyApp"), "."
            )),



            ###########
            # NotchKO #
            ###########
            tabItem(
                tabName = "dataset_notchKO",
                fluidRow(
                    style = "margin-top: -10px; margin-bottom: -10px; padding: 0px",
                    box(
                        title = "Celltype annotation",
                        plotOutput("annotation_plot", height = 290),
                        fluidRow(column(12,
                            align = "right", style = "margin-top: -295px; margin-right: -200",
                            downloadButton("download_annotation_plot", label = NULL)
                        )),
                        status = "primary", solidHeader = TRUE, width = 4
                    ),
                    box(
                        title = "Feature Expression",
                        plotOutput("feature_plot", height = 290),
                        fluidRow(column(12,
                            align = "right", style = "margin-top: -295px; margin-right: -200",
                            downloadButton("download_feature_plot", label = NULL)
                        )),
                        status = "primary", solidHeader = TRUE, width = 4
                    ),
                    box(
                        title = "Settings",
                        style = "height: 308px;",
                        div(
                            style = "margin-bottom: 0px",
                            selectizeInput("selected_gene", "Gene name",
                                choices = c(rownames(seurat), feature_mapping$fbid),
                                options = list(
                                    search_contains = TRUE,
                                    dropdownParent = "body",
                                    scoreThreshold = 1
                                )
                            )
                        ),
                        div(
                            style = "margin-top: 10px",
                            awesomeRadio("dim_reduction", HTML("Dimensionality reduction"), choices = c("UMAP", "PCA"), selected = "UMAP", inline = TRUE),
                            column(
                                style = "margin-left: -7px",
                                width = 3, offset = 0,
                                div(
                                    prettyCheckbox("order_cells", "Order cells", TRUE),
                                    prettyCheckbox("split_condition", "Split by condition", FALSE)
                                )
                            ),
                            column(
                                width = 1, offset = 7,
                                icon("info-circle", class = "icon-info", id = "order_cells_info"),
                                div(style = "margin-bottom: 15px;"),
                                icon("info-circle", class = "icon-info", id = "split_condition_info")
                            )
                        ),
                        actionButton("go", "Update All Plots!",
                            icon = icon("fa-solid fa-gears", class = "fa-lg"),
                            style = "font-weight: bold; background-color: #14E821; margin-top: 20px;"
                        ),
                        status = "info", solidHeader = TRUE, width = 3
                    ),
                    bsTooltip(id = "order_cells_info", title = "Should cells be plotted in order of expression?", placement = "right", trigger = "hover"),
                    bsTooltip(id = "split_condition_info", title = "Create seperate plot for Ctrl and Notch-knockout condition?", placement = "right", trigger = "hover"),
                ),
                div(
                    style = "margin-top: 10px",
                    conditionalPanel(
                        condition = "input.split_condition",
                        fluidRow(
                            box(
                                width = 4,
                                title = "Celltype annotation",
                                plotOutput("annotation_plot_split", height = 290),
                                fluidRow(
                                    column(12,
                                        align = "right", style = "margin-top: -295px; margin-right: -200",
                                        downloadButton("download_annotation_plot_split", label = NULL)
                                    )
                                ),
                                status = "primary", solidHeader = TRUE
                            ),
                            box(
                                width = 4, title = "Feature Expression",
                                plotOutput("feature_plot_split", height = 290),
                                fluidRow(
                                    column(12,
                                        align = "right", style = "margin-top: -295px; margin-right: -200",
                                        downloadButton("download_feature_plot_split", label = NULL)
                                    )
                                ),
                                status = "primary", solidHeader = TRUE
                            )
                        )
                    )
                ),
                fluidRow(
                    style = "margin-top: -10px; padding: 0px",
                    box(
                        title = "Expression across celltypes",
                        plotOutput("vln_expression", height = 200),
                        status = "primary", solidHeader = TRUE, width = 8,
                        fluidRow(column(12,
                            align = "right", style = "margin-top: -205px; margin-right: -200",
                            downloadButton("download_vln_expression", label = NULL)
                        )),
                    ),
                    box(
                        title = "Settings", status = "info", solidHeader = TRUE, width = 3,
                        div(
                            style = "column-count: 3; margin-bottom: 15px; ",
                            checkboxGroupInput("celltypes_vln_plot", "Available celltypes",
                                choices = c(
                                    "ISC", "EB", "EEP", "dEC", "daEC", "aEC",
                                    "mEC", "Copper", "LFC", "pEC", "EE", "MT"
                                ),
                                selected = c(
                                    "ISC", "EB", "EEP", "dEC", "aEC",
                                    "mEC", "Copper", "LFC", "pEC", "EE"
                                )
                            )
                        ),
                        div(
                            style = "margin-bottom: 0", # make box fit!
                            checkboxGroupInput("summary_stats", "Summary statistics",
                                choices = c("Include mean", "Include median"),
                                selected = "Include mean", inline = T
                            )
                        ),
                    )
                )
            ),
            ################
            # RNAi dataset #
            ################
            tabItem(
                tabName = "dataset_RNAi",
                fluidRow(
                    style = "margin-top: -10px; margin-bottom: -10px; padding: 0px",
                    box(
                        title = "Celltype annotation - Integrated Dataset",
                        plotOutput("annotation_plot_rnai", height = 290),
                        fluidRow(column(12,
                            align = "right", style = "margin-top: -295px; margin-right: -200",
                            downloadButton("download_annotation_plot_rnai", label = NULL)
                        )),
                        status = "primary", solidHeader = TRUE, width = 4
                    ),
                    box(
                        title = "Feature Expression - Integrated Dataset",
                        plotOutput("feature_plot_rnai", height = 290),
                        fluidRow(column(12,
                            align = "right", style = "margin-top: -295px; margin-right: -200",
                            downloadButton("download_feature_plot_rnai", label = NULL)
                        )),
                        status = "primary", solidHeader = TRUE, width = 4
                    ),
                    box(
                        title = "Settings",
                        div(
                            style = "margin-bottom: 0px; margin-top: 0px; padding: 0px;", # Set margin and padding to 0
                            selectizeInput("selected_gene_rnai", "Gene name",
                                choices = c(rownames(seurat_RNAi), feature_mapping$fbid),
                                options = list(
                                    search_contains = TRUE,
                                    dropdownParent = "body",
                                    scoreThreshold = 1
                                )
                            )
                        ),
                        awesomeRadio("dim_reduction_rnai", HTML("Dimensionality reduction"),
                            choices = c("UMAP", "PCA"), selected = "UMAP", inline = TRUE
                        ),
                        prettyCheckbox("order_cells_rnai", "Order cells", TRUE),
                        checkboxGroupInput("selected_conditions_rnai", "Select condition(s)",
                            choiceNames = list(
                                "Control",
                                HTML("<i>Notch<sup>RNAi</i></sup>"),
                                HTML("<i>Cph<sup>RNAi</i></sup>+<i>Notch<sup>RNAi</i></sup>"),
                                HTML("<i>Cph<sup>up</i></sup>")
                            ),
                            choiceValues = list("ctrl", "NotchRNAi", "NotchCphRNAi", "CphUp"),
                            inline = FALSE
                        ),
                        actionButton("go_rnai", "Update All Plots!",
                            icon = icon("fa-solid fa-gears", class = "fa-lg"),
                            style = "font-weight: bold; background-color: #14E821;"
                        ),
                        status = "info", solidHeader = TRUE, width = 3
                    )
                ),
                div(
                    conditionalPanel(
                        condition = "input.selected_conditions_rnai.includes('ctrl')",
                        fluidRow(
                            box(
                                width = 4,
                                title = HTML("Celltype annotation - Control"),
                                plotOutput("annotation_plot_ctrl_rnai", height = 290),
                                fluidRow(column(12,
                                    align = "right", style = "margin-top: -295px; margin-right: -200",
                                    downloadButton("download_annotation_plot_ctrl_rnai", label = NULL)
                                )),
                                status = "primary",
                                solidHeader = TRUE
                            ),
                            box(
                                width = 4,
                                title = HTML("Feature Expression - Control"),
                                plotOutput("feature_plot_ctrl_rnai", height = 290),
                                fluidRow(column(12,
                                    align = "right", style = "margin-top: -295px; margin-right: -200",
                                    downloadButton("download_feature_plot_ctrl_rnai", label = NULL)
                                )),
                                status = "primary",
                                solidHeader = TRUE
                            )
                        )
                    )
                ),
                div(
                    conditionalPanel(
                        condition = "input.selected_conditions_rnai.includes('NotchRNAi')",
                        fluidRow(
                            box(
                                width = 4,
                                title = HTML("Celltype annotation - <i>Notch<sup>RNAi</i></sup>"),
                                plotOutput("annotation_plot_NotchRNAi_rnai", height = 290),
                                fluidRow(column(12,
                                    align = "right", style = "margin-top: -295px; margin-right: -200",
                                    downloadButton("download_annotation_plot_NotchRNAi_rnai", label = NULL)
                                )),
                                status = "primary",
                                solidHeader = TRUE
                            ),
                            box(
                                width = 4,
                                title = HTML("Feature Expression - <i>Notch<sup>RNAi</i></sup>"),
                                plotOutput("feature_plot_NotchRNAi_rnai", height = 290),
                                fluidRow(column(12,
                                    align = "right", style = "margin-top: -295px; margin-right: -200",
                                    downloadButton("download_feature_plot_NotchRNAi_rnai", label = NULL)
                                )),
                                status = "primary",
                                solidHeader = TRUE
                            )
                        )
                    )
                ),
                div(
                    conditionalPanel(
                        condition = "input.selected_conditions_rnai.includes('NotchCphRNAi')",
                        fluidRow(
                            box(
                                width = 4,
                                title = HTML("Celltype annotation - <i>Cph<sup>RNAi</i></sup>+<i>Notch<sup>RNAi</i></sup>"),
                                plotOutput("annotation_plot_NotchCphRNAi_rnai", height = 290),
                                fluidRow(column(12,
                                    align = "right", style = "margin-top: -295px; margin-right: -200",
                                    downloadButton("download_annotation_plot_NotchCphRNAi_rnai", label = NULL)
                                )),
                                status = "primary",
                                solidHeader = TRUE
                            ),
                            box(
                                width = 4,
                                title = HTML("Feature Expression - <i>Cph<sup>RNAi</i></sup>+<i>Notch<sup>RNAi</i></sup>"),
                                plotOutput("feature_plot_NotchCphRNAi_rnai", height = 290),
                                fluidRow(column(12,
                                    align = "right", style = "margin-top: -295px; margin-right: -200",
                                    downloadButton("download_feature_plot_NotchCphRNAi_rnai", label = NULL)
                                )),
                                status = "primary",
                                solidHeader = TRUE
                            )
                        )
                    )
                ),
                div(
                    conditionalPanel(
                        condition = "input.selected_conditions_rnai.includes('CphUp')",
                        fluidRow(
                            box(
                                width = 4,
                                title = HTML("Celltype annotation - <i>Cph<sup>up</i></sup>"),
                                plotOutput("annotation_plot_CphUp_rnai", height = 290),
                                fluidRow(column(12,
                                    align = "right", style = "margin-top: -295px; margin-right: -200",
                                    downloadButton("download_annotation_plot_CphUp_rnai", label = NULL)
                                )),
                                status = "primary",
                                solidHeader = TRUE
                            ),
                            box(
                                width = 4,
                                title = HTML("Feature Expression - <i>Cph<sup>up</i></sup>"),
                                plotOutput("feature_plot_CphUp_rnai", height = 290),
                                fluidRow(column(12,
                                    align = "right", style = "margin-top: -295px; margin-right: -200",
                                    downloadButton("download_feature_plot_CphUp_rnai", label = NULL)
                                )),
                                status = "primary",
                                solidHeader = TRUE
                            )
                        )
                    )
                ),
                fluidRow(
                    style = "margin-top: -10px; padding: 0px",
                    box(
                        title = "Expression across celltypes",
                        plotOutput("vln_expression_rnai", height = 200),
                        fluidRow(column(12,
                            align = "right", style = "margin-top: -205px; margin-right: -200",
                            downloadButton("download_vln_expression_rnai", label = NULL)
                        )),
                        status = "primary",
                        solidHeader = TRUE, width = 8,
                    ),
                    box(
                        title = "Settings", solidHeader = TRUE, width = 3,
                        status = "info",
                        div(
                            style = "column-count: 3; margin-bottom: 15px",
                            checkboxGroupInput("celltypes_vln_plot_rnai", "Available celltypes",
                                choices = c(
                                    "ISC", "EB", "EEP", "dEC", "daEC", "aEC",
                                    "mEC", "Copper", "LFC", "pEC", "EE", "MT"
                                ),
                                selected = c(
                                    "ISC", "EB", "EEP", "dEC", "aEC",
                                    "mEC", "Copper", "LFC", "pEC", "EE"
                                )
                            )
                        ),
                        div(
                            style = "margin-bottom: 0px;",
                            checkboxGroupInput("summary_stats_rnai", "Summary statistics",
                                choices = c("Include mean", "Include median"),
                                selected = "Include mean", inline = T
                            )
                        )
                    )
                )
            ),
            #######################
            # REGIONAL EXPRESSION #
            #######################
            tabItem(
                tabName = "regional_expression",
                fluidRow(
                    style = "margin-top: -10px; margin-bottom: -10px; padding: 0px",
                    box(
                        title = "Gene Expression across Midgut Regions",
                        plotOutput("regional_plot", height = 400),
                        fluidRow(column(12,
                            align = "right", style = "margin-top: -405px; margin-right: -200",
                            downloadButton("download_regional_plot", label = NULL)
                        )),
                        status = "primary", solidHeader = TRUE, width = 8
                    ),
                    box(
                        title = "Settings",
                        style = "height: 418px;",
                        div(
                            style = "margin-bottom: 0px",
                            selectizeInput("selected_gene_regional", "Gene name",
                                choices = c(rownames(seurat), feature_mapping$fbid),
                                options = list(
                                    search_contains = TRUE,
                                    dropdownParent = "body",
                                    scoreThreshold = 1
                                )
                            )
                        ),
                        actionButton("go_regional", "Update Plot!",
                            icon = icon("fa-solid fa-gears", class = "fa-lg"),
                            style = "font-weight: bold; background-color: #14E821; margin-top: 20px;"
                        ),
                        status = "info", solidHeader = TRUE, width = 4
                    )
                )
            ),
            ###############
            # CONTACT TAB #
            ###############
            tabItem(
                tabName = "contact",
                h2("Contact Information"),
                p("If you have any questions about the publication or Shiny app, feel free to reach out to any of the people listed below \U1F680."),
                p(
                    "In case you find any bugs or have feature requests please post them on ",
                    a("GitHub", href = "https://github.com/nickhir/Chronophage/issues"),
                    "."
                ),
                div(style = "margin-bottom: 20px;"), # Add a margin to create space
                h4("Joint first authors", style = "font-weight: bold;"),
                p("Siamak Redhai,",
                    a("siamak.redhai@dkfz-heidelberg.de", href = "mailto:siamak.redhai@dkfz-heidelberg.de"),
                    ",",
                    a("German Cancer Research Center (DKFZ)", href = "https://www.dkfz.de/en/index.html"),
                    style = "margin-bottom: 3px;"
                ),
                p("Nick Hirschmüller,",
                    a("hirschmueller.nick@gmail.com", href = "mailto:hirschmueller.nick@gmail.com"),
                    ",",
                    a("University of Cambridge", href = "https://www.cam.ac.uk/"),
                    style = "margin-bottom: 3px;"
                ),
                h4("Corresponding authors", style = "font-weight: bold;"),
                p("Siamak Redhai,",
                    a("siamak.redhai@dkfz-heidelberg.de", href = "mailto:siamak.redhai@dkfz-heidelberg.de"),
                    ",",
                    a("German Cancer Research Center (DKFZ)", href = "https://www.dkfz.de/en/index.html"),
                    style = "margin-bottom: 3px;"
                ),
                p("Michael Boutros,",
                    a("m.boutros@dkfz-heidelberg.de", href = "mailto:m.boutros@dkfz-heidelberg.de"),
                    ",",
                    a("German Cancer Research Center (DKFZ)", href = "https://www.dkfz.de/en/index.html"),
                    style = "margin-bottom: 3px;"
                ),
                p("Wolfgang Huber,",
                    a("wolfgang.huber@embl.org", href = "mailto:wolfgang.huber@embl.org"),
                    ",",
                    a("European Molecular Biology Laboratory (EMBL)", href = "https://www.embl.org/"),
                    style = "margin-bottom: 3px;"
                ),
            )
        )
    )
)





##########
# SERVER #
##########
server <- function(input, output, session) {
    # server side selection
    updateSelectizeInput(session, "selected_gene",
        choices = c(rownames(seurat), feature_mapping$fbid), server = TRUE,
        selected = "Cph"
    )

    updateSelectizeInput(session, "selected_gene_rnai",
        choices = c(rownames(seurat), feature_mapping$fbid), server = TRUE,
        selected = "esg"
    )

    updateSelectizeInput(session, "selected_gene_regional",
        choices = c(rownames(seurat), feature_mapping$fbid), server = TRUE,
        selected = "esg"
    )


    # Helper function to create initial plots
    create_default_plots <- function(srt, gene, dataset_type = "notch") {
        default_celltypes <- c("ISC", "EB", "EEP", "dEC", "daEC", "aEC", "mEC", "Copper", "LFC", "pEC", "EE")

        scale_transform <- if (dataset_type == "rnai") "scale_x_reverse" else "scale_y_reverse"

        list(
            annotation = annotation_plot_fun(srt, "UMAP", "combined") %>% apply_scale_transform(scale_transform),
            feature = feature_plot_fun(srt, gene, "UMAP", TRUE, "combined") %>% apply_scale_transform(scale_transform),
            violin = jitter_plot_fun(srt, gene, default_celltypes, "Include mean", "combined")
        )
    }

    # By default, we load some standard plots. So the user sees something upon startup
    notch_plots <- create_default_plots(seurat, "Cph", "notch")
    output$annotation_plot <- renderPlot(notch_plots$annotation)
    output$feature_plot <- renderPlot(notch_plots$feature)
    output$vln_expression <- renderPlot(notch_plots$violin)

    rnai_plots <- create_default_plots(seurat_RNAi, "esg", "rnai")
    output$annotation_plot_rnai <- renderPlot(rnai_plots$annotation)
    output$feature_plot_rnai <- renderPlot(rnai_plots$feature)
    output$vln_expression_rnai <- renderPlot(rnai_plots$violin)

    output$regional_plot <- renderPlot(vln_gene_by_region(seurat, "esg"))

    ##################
    # NOTCH KO PLOTS #
    ##################
    observeEvent(input$go, {
        dim_red <- input$dim_reduction
        order_cells <- input$order_cells
        celltypes <- input$celltypes_vln_plot
        gene_name <- map_gene_name(input$selected_gene, feature_mapping)
        summary_stats <- input$summary_stats
        split_condition <- input$split_condition

        if (!split_condition) {
            # Combined plots
            output$annotation_plot <- renderPlot(
                annotation_plot_fun(seurat, dim_red, "combined") %>%
                    apply_scale_transform("scale_y_reverse")
            ) %>% bindCache(dim_red, "combined", cache = "session")

            output$feature_plot <- renderPlot(
                feature_plot_fun(seurat, gene_name, dim_red, order_cells, "combined") %>%
                    apply_scale_transform("scale_y_reverse")
            ) %>% bindCache(dim_red, gene_name, order_cells, "combined", cache = "session")

            output$vln_expression <- renderPlot(
                jitter_plot_fun(seurat, gene_name, celltypes, summary_stats, "combined")
            ) %>% bindCache(celltypes, gene_name, summary_stats, "combined", cache = "session")
        }

        # Split condition plots
        if (split_condition) {
            output$annotation_plot <- renderPlot(
                annotation_plot_fun(seurat, dim_red, "ctrl") %>%
                    apply_scale_transform("scale_y_reverse")
            ) %>% bindCache(dim_red, "ctrl", cache = "session")

            output$annotation_plot_split <- renderPlot(
                annotation_plot_fun(seurat, dim_red, "notch") %>%
                    apply_scale_transform("scale_y_reverse")
            ) %>% bindCache(dim_red, "notch", cache = "session")

            output$feature_plot <- renderPlot(
                feature_plot_fun(seurat, gene_name, dim_red, order_cells, "ctrl") %>%
                    apply_scale_transform("scale_y_reverse")
            ) %>% bindCache(dim_red, gene_name, order_cells, "ctrl", cache = "session")

            output$feature_plot_split <- renderPlot(
                feature_plot_fun(seurat, gene_name, dim_red, order_cells, "notch") %>%
                    apply_scale_transform("scale_y_reverse")
            ) %>% bindCache(dim_red, gene_name, order_cells, "notch", cache = "session")

            output$vln_expression <- renderPlot(
                jitter_plot_fun(seurat, gene_name, celltypes, summary_stats, "split")
            ) %>% bindCache(celltypes, gene_name, summary_stats, "split", cache = "session")
        }
    })

    ##############
    # RNAi PLOTS #
    ##############
    observeEvent(input$go_rnai, {
        dim_red_rnai <- input$dim_reduction_rnai
        order_cells_rnai <- input$order_cells_rnai
        celltypes_rnai <- input$celltypes_vln_plot_rnai
        gene_name_rnai <- map_gene_name(input$selected_gene_rnai, feature_mapping)
        summary_stats_rnai <- input$summary_stats_rnai
        selected_conditions_rnai <- input$selected_conditions_rnai

        # Combined plots
        output$annotation_plot_rnai <- renderPlot(
            annotation_plot_fun(seurat_RNAi, dim_red_rnai, "combined") %>%
                apply_scale_transform("scale_x_reverse")
        ) %>% bindCache(dim_red_rnai, "combined", cache = "session")

        output$feature_plot_rnai <- renderPlot(
            feature_plot_fun(seurat_RNAi, gene_name_rnai, dim_red_rnai, order_cells_rnai, "combined") %>%
                apply_scale_transform("scale_x_reverse")
        ) %>% bindCache(dim_red_rnai, gene_name_rnai, order_cells_rnai, "combined", cache = "session")

        output$vln_expression_rnai <- renderPlot(
            jitter_plot_fun(seurat_RNAi, gene_name_rnai, celltypes_rnai, summary_stats_rnai, "combined")
        ) %>% bindCache(celltypes_rnai, gene_name_rnai, summary_stats_rnai, "combined", cache = "session")

        # Condition-specific plots
        if (!is.null(selected_conditions_rnai)) {
            # Helper function to create condition-specific plots
            create_condition_plots <- function(condition, output_suffix) {
                output[[paste0("annotation_plot_", condition, "_rnai")]] <- renderPlot(
                    annotation_plot_fun(seurat_RNAi, dim_red_rnai, condition) %>%
                        apply_scale_transform("scale_x_reverse")
                ) %>% bindCache(dim_red_rnai, condition, cache = "session")

                output[[paste0("feature_plot_", condition, "_rnai")]] <- renderPlot(
                    feature_plot_fun(seurat_RNAi, gene_name_rnai, dim_red_rnai, order_cells_rnai, condition) %>%
                        apply_scale_transform("scale_x_reverse")
                ) %>% bindCache(dim_red_rnai, gene_name_rnai, order_cells_rnai, condition, cache = "session")
            }

            # Create plots for each selected condition
            conditions_map <- list(
                "ctrl" = "ctrl",
                "NotchRNAi" = "NotchRNAi",
                "CphUp" = "CphUp",
                "NotchCphRNAi" = "NotchCphRNAi"
            )

            for (condition in selected_conditions_rnai) {
                if (condition %in% names(conditions_map)) {
                    create_condition_plots(condition, conditions_map[[condition]])
                }
            }

            # Update violin plot with selected conditions
            output$vln_expression_rnai <- renderPlot(
                jitter_plot_fun_rnai(seurat_RNAi, gene_name_rnai, celltypes_rnai, summary_stats_rnai, selected_conditions_rnai)
            ) %>% bindCache(celltypes_rnai, gene_name_rnai, summary_stats_rnai, selected_conditions_rnai, cache = "session")
        }
    })

    #########################
    # REGIONAL EXPRESSION #
    #########################
    observeEvent(input$go_regional, {
        gene_name_regional <- map_gene_name(input$selected_gene_regional, feature_mapping)

        output$regional_plot <- renderPlot(
            vln_gene_by_region(seurat, gene_name_regional)
        ) %>% bindCache(gene_name_regional, cache = "session")
    })

    #################
    # DownloadCalls #
    #################
    # Helper function to create download handlers
    create_download_handler <- function(plot_func, filename_prefix, file_ext = ".pdf", use_ggsave = FALSE, use_rnai_inputs = FALSE) {
        downloadHandler(
            filename = function() {
                # Use context-aware input selection based on use_rnai_inputs flag
                if (use_rnai_inputs) {
                    gene <- if (exists("input") && !is.null(input$selected_gene_rnai)) {
                        map_gene_name(input$selected_gene_rnai, feature_mapping)
                    } else {
                        "gene"
                    }
                    dim_red <- if (exists("input") && !is.null(input$dim_reduction_rnai)) {
                        input$dim_reduction_rnai
                    } else {
                        "UMAP"
                    }
                } else {
                    gene <- if (exists("input") && !is.null(input$selected_gene)) {
                        map_gene_name(input$selected_gene, feature_mapping)
                    } else {
                        "gene"
                    }
                    dim_red <- if (exists("input") && !is.null(input$dim_reduction)) {
                        input$dim_reduction
                    } else {
                        "UMAP"
                    }
                }

                paste0(filename_prefix, "_", if (grepl("annotation", filename_prefix)) dim_red else gene, file_ext)
            },
            content = function(file) {
                p <- plot_func()
                if (use_ggsave) {
                    ggsave(filename = file, width = 7, height = 2.5, plot = p, scale = 1.4)
                } else {
                    pdf(file, width = 7.31, height = 5.61)
                    print(p)
                    dev.off()
                }
            }
        )
    }

    observe({
        gene_name <- map_gene_name(input$selected_gene, feature_mapping)
        split_condition <- input$split_condition

        if (!split_condition) {
            # Combined condition downloads
            output$download_annotation_plot <- create_download_handler(
                function() annotation_plot_fun(seurat, input$dim_reduction, "combined"),
                "annotation"
            )

            output$download_feature_plot <- create_download_handler(
                function() feature_plot_fun(seurat, gene_name, input$dim_reduction, input$order_cells, "combined"),
                "feature_expression"
            )

            output$download_vln_expression <- create_download_handler(
                function() jitter_plot_fun(seurat, gene_name, input$celltypes_vln_plot, input$summary_stats, "combined"),
                "celltype_expression", use_ggsave = TRUE
            )
        } else {
            # Split condition downloads
            output$download_annotation_plot <- create_download_handler(
                function() annotation_plot_fun(seurat, input$dim_reduction, "ctrl"),
                "annotation"
            )

            output$download_feature_plot <- create_download_handler(
                function() feature_plot_fun(seurat, gene_name, input$dim_reduction, input$order_cells, "ctrl"),
                "feature_expression"
            )

            output$download_vln_expression <- create_download_handler(
                function() jitter_plot_fun(seurat, gene_name, input$celltypes_vln_plot, input$summary_stats, "split"),
                "celltype_expression", use_ggsave = TRUE
            )

            output$download_annotation_plot_split <- create_download_handler(
                function() annotation_plot_fun(seurat, input$dim_reduction, "notch"),
                "annotation_notch"
            )

            output$download_feature_plot_split <- create_download_handler(
                function() feature_plot_fun(seurat, gene_name, input$dim_reduction, input$order_cells, "notch"),
                "feature_expression_notch"
            )
        }
    })

    # RNAi download handlers
    observe({
        gene_name_rnai <- map_gene_name(input$selected_gene_rnai, feature_mapping)

        # Combined RNAi dataset downloads
        output$download_annotation_plot_rnai <- create_download_handler(
            function() annotation_plot_fun(seurat_RNAi, input$dim_reduction_rnai, "combined") %>%
                apply_scale_transform("scale_x_reverse"),
            "annotation_rnai", use_rnai_inputs = TRUE
        )

        output$download_feature_plot_rnai <- create_download_handler(
            function() feature_plot_fun(seurat_RNAi, gene_name_rnai, input$dim_reduction_rnai, input$order_cells_rnai, "combined") %>%
                apply_scale_transform("scale_x_reverse"),
            "feature_expression_rnai", use_rnai_inputs = TRUE
        )

        output$download_vln_expression_rnai <- create_download_handler(
            function() {
                if (!is.null(input$selected_conditions_rnai)) {
                    jitter_plot_fun_rnai(seurat_RNAi, gene_name_rnai, input$celltypes_vln_plot_rnai, input$summary_stats_rnai, input$selected_conditions_rnai)
                } else {
                    jitter_plot_fun(seurat_RNAi, gene_name_rnai, input$celltypes_vln_plot_rnai, input$summary_stats_rnai, "combined")
                }
            },
            "celltype_expression_rnai", use_ggsave = TRUE, use_rnai_inputs = TRUE
        )

        # Condition-specific downloads
        # Control condition
        output$download_annotation_plot_ctrl_rnai <- create_download_handler(
            function() annotation_plot_fun(seurat_RNAi, input$dim_reduction_rnai, "ctrl") %>%
                apply_scale_transform("scale_x_reverse"),
            "annotation_ctrl_rnai", use_rnai_inputs = TRUE
        )

        output$download_feature_plot_ctrl_rnai <- create_download_handler(
            function() feature_plot_fun(seurat_RNAi, gene_name_rnai, input$dim_reduction_rnai, input$order_cells_rnai, "ctrl") %>%
                apply_scale_transform("scale_x_reverse"),
            "feature_expression_ctrl_rnai", use_rnai_inputs = TRUE
        )

        # NotchRNAi condition
        output$download_annotation_plot_NotchRNAi_rnai <- create_download_handler(
            function() annotation_plot_fun(seurat_RNAi, input$dim_reduction_rnai, "NotchRNAi") %>%
                apply_scale_transform("scale_x_reverse"),
            "annotation_NotchRNAi_rnai", use_rnai_inputs = TRUE
        )

        output$download_feature_plot_NotchRNAi_rnai <- create_download_handler(
            function() feature_plot_fun(seurat_RNAi, gene_name_rnai, input$dim_reduction_rnai, input$order_cells_rnai, "NotchRNAi") %>%
                apply_scale_transform("scale_x_reverse"),
            "feature_expression_NotchRNAi_rnai", use_rnai_inputs = TRUE
        )

        # NotchCphRNAi condition
        output$download_annotation_plot_NotchCphRNAi_rnai <- create_download_handler(
            function() annotation_plot_fun(seurat_RNAi, input$dim_reduction_rnai, "NotchCphRNAi") %>%
                apply_scale_transform("scale_x_reverse"),
            "annotation_NotchCphRNAi_rnai", use_rnai_inputs = TRUE
        )

        output$download_feature_plot_NotchCphRNAi_rnai <- create_download_handler(
            function() feature_plot_fun(seurat_RNAi, gene_name_rnai, input$dim_reduction_rnai, input$order_cells_rnai, "NotchCphRNAi") %>%
                apply_scale_transform("scale_x_reverse"),
            "feature_expression_NotchCphRNAi_rnai", use_rnai_inputs = TRUE
        )

        # CphUp condition
        output$download_annotation_plot_CphUp_rnai <- create_download_handler(
            function() annotation_plot_fun(seurat_RNAi, input$dim_reduction_rnai, "CphUp") %>%
                apply_scale_transform("scale_x_reverse"),
            "annotation_CphUp_rnai", use_rnai_inputs = TRUE
        )

        output$download_feature_plot_CphUp_rnai <- create_download_handler(
            function() feature_plot_fun(seurat_RNAi, gene_name_rnai, input$dim_reduction_rnai, input$order_cells_rnai, "CphUp") %>%
                apply_scale_transform("scale_x_reverse"),
            "feature_expression_CphUp_rnai", use_rnai_inputs = TRUE
        )
    })

    # Regional expression download handler
    observe({
        gene_name_regional <- map_gene_name(input$selected_gene_regional, feature_mapping)

        output$download_regional_plot <- create_download_handler(
            function() vln_gene_by_region(seurat, gene_name_regional),
            "regional_expression", use_ggsave = TRUE
        )
    })
}

shinyApp(ui, server)
