library(shiny)
library(shinydashboard)
library(magrittr)
library(ggplot2)
library(shinyWidgets)
library(shinyBS)
library(dplyr)
library(Seurat)
library(ggtext)

options(styler.addins_style_transformer = "styler::tidyverse_style(indent_by = 4L)")

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

    # generate the plot
    plot_data <- srt@meta.data[, c("celltype_manual", dimensions, "perturbation")] %>%
        dplyr::sample_frac()

    if (split_conditions == "combined") {
        p <- annotation_plot_call(plot_data, dimensions) +
            ggtitle("Combined Dataset")
        if (n_distinct(srt$perturbation) == 2) {
            p <- p + scale_x_reverse()
        }
    } else if (split_conditions == "ctrl") {
        p <- annotation_plot_call(plot_data %>% filter(perturbation == "ctrl"), dimensions) +
            ggtitle("Control Condition") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))

        if (n_distinct(srt$perturbation) == 2) {
            p <- p + scale_x_reverse()
        }
    } else if (split_conditions == "notch") {
        p <- annotation_plot_call(plot_data %>% filter(perturbation == "notch"), dimensions) +
            ggtitle("*Notch*<sup><i>sgRNAx2</i></sup>") +
            theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15)) +
            scale_x_reverse()
    } else if (split_conditions == "NotchRNAi") {
        p <- annotation_plot_call(plot_data %>% filter(perturbation == "NotchRNAi"), dimensions) +
            ggtitle("*Notch*<sup><i>RNAi</i></sup>") +
            theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
    } else if (split_conditions == "NotchCphRNAi") {
        p <- annotation_plot_call(plot_data %>% filter(perturbation == "NotchCphRNAi"), dimensions) +
            ggtitle("*Cph*<sup><i>RNAi</i></sup>+*Notch*<sup><i>RNAi</i></sup>") +
            theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
    } else if (split_conditions == "CphUp") {
        p <- annotation_plot_call(plot_data %>% filter(perturbation == "CphUp"), dimensions) +
            ggtitle("*Cph*<sup>up</sup>") +
            theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
    }
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
    mdata <- srt@meta.data[, c("celltype_manual", dimensions, "perturbation")]
    exp_data <- FetchData(srt, gene)
    colnames(exp_data) <- gene
    stopifnot(all(rownames(mdata) == rownames(exp_data)))

    plot_data <- cbind(mdata, exp_data)

    if (order_cells) {
        plot_data <- plot_data %>%
            dplyr::arrange(., !!sym(gene))
    } else {
        plot_data <- plot_data %>%
            dplyr::sample_frac()
    }

    if (split_conditions == "combined") {
        p <- feature_plot_call(plot_data, dimensions, gene) +
            ggtitle("Combined Dataset")
        if (n_distinct(srt$perturbation) == 2) {
            p <- p + scale_x_reverse()
        }
    } else if (split_conditions == "ctrl") {
        p <- feature_plot_call(plot_data %>% filter(perturbation == "ctrl"), dimensions, gene) +
            ggtitle("Control Condition") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
        if (n_distinct(srt$perturbation) == 2) {
            p <- p + scale_x_reverse()
        }
    } else if (split_conditions == "notch") {
        p <- feature_plot_call(plot_data %>% filter(perturbation == "notch"), dimensions, gene) +
            ggtitle("*Notch*<sup><i>sgRNAx2</i></sup>") +
            theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15)) +
            scale_x_reverse()
    } else if (split_conditions == "NotchRNAi") {
        p <- feature_plot_call(plot_data %>% filter(perturbation == "NotchRNAi"), dimensions, gene) +
            ggtitle("*Notch*<sup><i>RNAi</i></sup>") +
            theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
    } else if (split_conditions == "NotchCphRNAi") {
        p <- feature_plot_call(plot_data %>% filter(perturbation == "NotchCphRNAi"), dimensions, gene) +
            ggtitle("*Cph*<sup><i>RNAi</i></sup>+*Notch*<sup><i>RNAi</i></sup>") +
            theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
    } else if (split_conditions == "CphUp") {
        p <- feature_plot_call(plot_data %>% filter(perturbation == "NotchCphRNAi"), dimensions, gene) +
            ggtitle("*Cph*<sup>Up</sup>") +
            theme(plot.title = element_markdown(hjust = 0.5, face = "bold", size = 15))
    }
    return(p)
}

jitter_plot_fun <- function(srt, gene, celltypes, summary_stats, split_conditions) {
  mdata <- srt@meta.data[, c("celltype_manual", "perturbation")]
  exp_data <- FetchData(srt, gene)
  colnames(exp_data) <- gene
  stopifnot(all(rownames(mdata) == rownames(exp_data)))
  
  plot_data <- cbind(mdata, exp_data) %>%
    filter(celltype_manual %in% celltypes)
  
  if (split_conditions == "combined") {
    p <- ggplot(plot_data, aes(x = celltype_manual, y = !!sym(gene), color = celltype_manual, fill = celltype_manual)) +
      geom_jitter(size = 0.65, stroke = 0.6, width = 0.2, shape = 21) +
      geom_violin(
        adjust = 1, trim = TRUE, color = "black", show.legend = F,
        scale = "width", fill = "#00000000",
        width = 0.72 # adjust width so the violine does not stick out too much
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
    
    if (!is.null(summary_stats) & "Include mean" %in% summary_stats) {
      p <- p +
        stat_summary(fun = "mean", geom = "point", color = "black", size = 3, shape = 8, stroke = 1.2)
    }
    if (!is.null(summary_stats) & "Include median" %in% summary_stats) {
      p <- p +
        stat_summary(fun = "median", geom = "point", color = "black", size = 3, shape = 2, stroke = 1.2)
    }
  } else if (split_conditions == "split") {
    p <- ggplot(plot_data, aes(x = celltype_manual, y = !!sym(gene))) +
      geom_point(
        mapping = aes(color = perturbation, fill = perturbation), size = 0.65, stroke = 0.4, shape = 21,
        position = position_jitterdodge()
      ) +
      geom_violin(
        mapping = aes(fill = perturbation),
        adjust = 1, trim = TRUE,
        scale = "width", show.legend = F,
        width = 0.75 # adjust width so the violine does not stick out too much
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
    
    if (!is.null(summary_stats) & "Include mean" %in% summary_stats) {
      p <- p +
        stat_summary(
          fun = "mean", geom = "point", aes(group = perturbation),
          color = "black", size = 3, shape = 8, stroke = 1.2,
          position = position_dodge2(width = 0.75)
        )
    }
    if (!is.null(summary_stats) & "Include median" %in% summary_stats) {
      p <- p +
        stat_summary(
          fun = "median", geom = "point", aes(group = perturbation),
          color = "black", size = 3, shape = 2, stroke = 1.2,
          position = position_dodge2(width = 0.75)
        )
    }
    if (n_distinct(srt$perturbation) == 3) {
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
  
  
  return(p +
           ggtitle(gene))
}


jitter_plot_fun_rnai <- function(srt, gene, celltypes, summary_stats, conditions) {
  mdata <- srt@meta.data[, c("celltype_manual", "perturbation")]
  exp_data <- FetchData(srt, gene)
  colnames(exp_data) <- gene
  stopifnot(all(rownames(mdata) == rownames(exp_data)))
  
  plot_data <- cbind(mdata, exp_data) %>%
    filter(celltype_manual %in% celltypes) %>% 
    filter(perturbation %in% conditions)
  
  p <- ggplot(plot_data, aes(x = celltype_manual, y = !!sym(gene))) +
    geom_point(
      mapping = aes(color = perturbation, fill = perturbation), size = 0.65, stroke = 0.4, shape = 21,
      position = position_jitterdodge()
    ) +
    geom_violin(
      mapping = aes(fill = perturbation),
      adjust = 1, trim = TRUE,
      scale = "width", show.legend = F,
      width = 0.75 # adjust width so the violine does not stick out too much
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
  
  if (!is.null(summary_stats) & "Include mean" %in% summary_stats) {
    p <- p +
      stat_summary(
        fun = "mean", geom = "point", aes(group = perturbation),
        color = "black", size = 3, shape = 8, stroke = 1.2,
        position = position_dodge2(width = 0.75)
      )
  }
  if (!is.null(summary_stats) & "Include median" %in% summary_stats) {
    p <- p +
      stat_summary(
        fun = "median", geom = "point", aes(group = perturbation),
        color = "black", size = 3, shape = 2, stroke = 1.2,
        position = position_dodge2(width = 0.75)
      )
  }
  p <- p + 
    scale_color_manual(values=color_mapping_rnai)+
    scale_fill_manual(values = alpha(color_mapping_rnai, 0))
  
  return(p +
           ggtitle(gene))
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
            #############
            # ABOUT TAB #
            #############
            tabItem(
                tabName = "about",
                h4("Intestinal stem cell proliferation and differentiation depends on the self-repressive zinc finger transcription factor BCL11/Cph",
                    align = "center",
                    style = "width: 80%;  margin: 0 auto; font-weight: bold; font-style: italic;"
                ),
                div(style = "margin-bottom: 10px;"), # Add a margin to create space

                h5(HTML("Siamak Redhai*<sup>#</sup>, Nick Hirschmüller*, Tianyu Wang, Shivohum Bahuguna,
              Svenja Leible, Stefan Peidli, Erica Valentani, Sviat Kharuk, Michaela Holemz, Lea Bräckow,
                        Fillip Port, David Ibberson, Wolfgang Huber<sup>#</sup>, Michael Boutros<sup>#</sup>"),
                    align = "center",
                    style = "width: 60%;  margin: 0 auto; color: darkgrey;"
                ),
                div(style = "margin-bottom: 5px;"), # Add a margin to create space

                h6("*Contributed equally | #Co-Corresponding",
                    align = "center",
                    style = "width: 80%;  margin: 0 auto; color: darkgrey;"
                ),
                div(style = "margin-bottom: 30px;"), # Add a margin to create space
                h3("Abstract", style = "font-weight: bold"),
                p("The molecular programs that drive proliferation and differentiation of intestinal stem cells (ISCs) are essential for organismal fitness. Notch signalling regulates the binary fate decision of ISCs, favouring enterocyte commitment when Notch activity is high and enteroendocrine cell (EE) fate when activity is low. However, the temporal dynamics of transcription factor expression and its involvement in generating distinct EE lineages across different intestinal regions remain largely uncharacterised. Here, we find that the expression of the C2H2-type zinc-finger transcription factor Chronophage (Cph), homologous to mammalian BCL11, increases early along the ISC-to-EE lineage when Notch is inactivated. We show that the expression of Cph is regulated by the Achaete-Scute Complex (AS-C) gene, scute, which directly binds to multiple sites within the Cph locus to promote its expression. Our genetic and single-cell RNA sequencing experiments demonstrate that Cph maintains the ISC and EE populations across different regions and is required to remodel the transcriptome of ISCs with low Notch activity. By performing in vivo chromatin profiling, we find that Cph activates the expression of key genes involved in proliferation and lineage commitment, while simultaneously repressing its own expression. This inhibitory feedback loop prevents ISCs from undergoing autophagy-dependent cell death and ensures differentiation is faithfully executed. Our findings shed light on the mechanism by which Cph sustains intestinal epithelial homeostasis and could represent a conserved strategy for balancing proliferation and differentiation in other tissues and species."),
                p(
                    "This Shiny App allows users to interactively explore the datasets accompanying the publication.",
                    "The source code for the app can be found on ", a("GitHub", href = "https://github.com/nickhir/IntestiMap"), "."
                )
            ),



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
                        status = "primary", solidHeader = TRUE, width = 4
                    ),
                    box(
                        title = "Feature Expression - Integrated Dataset",
                        plotOutput("feature_plot_rnai", height = 290),
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
                                status = "primary",
                                solidHeader = TRUE
                            ),
                            box(
                                width = 4,
                                title = HTML("Feature Expression - Control"),
                                plotOutput("feature_plot_ctrl_rnai", height = 290),
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
                                status = "primary",
                                solidHeader = TRUE
                            ),
                            box(
                                width = 4,
                                title = HTML("Feature Expression - <i>Notch<sup>RNAi</i></sup>"),
                                plotOutput("feature_plot_NotchRNAi_rnai", height = 290),
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
                                status = "primary",
                                solidHeader = TRUE
                            ),
                            box(
                                width = 4,
                                title = HTML("Feature Expression - <i>Cph<sup>RNAi</i></sup>+<i>Notch<sup>RNAi</i></sup>"),
                                plotOutput("feature_plot_NotchCphRNAi_rnai", height = 290),
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
                                status = "primary",
                                solidHeader = TRUE
                            ),
                            box(
                                width = 4,
                                title = HTML("Feature Expression - <i>Cph<sup>up</i></sup>"),
                                plotOutput("feature_plot_CphUp_rnai", height = 290),
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
            ###############
            # CONTACT TAB #
            ###############
            tabItem(
                tabName = "contact",
                h2("Contact Information"),
                p("If you have any questions about the publication or Shiny app, feel free to reach out to any of the people listed below \U1F680."),
                p(
                    "In case you find any bugs or have feature requests please post them on ",
                    a("GitHub", href = "https://github.com/nickhir/IntestiMap/issues"),
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
                    a("nh608@cam.ac.uk", href = "mailto:nh608@cam.ac.uk"),
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


    # By default, we load some standard plots. So the user sees something upon startup
    output$annotation_plot <- renderPlot(annotation_plot_fun(seurat,
        dim_red = "UMAP",
        split_conditions = "combined"
    ) +
        ggtitle("Combined Dataset")+
      scale_y_reverse())

    output$feature_plot <- renderPlot(feature_plot_fun(seurat,
        gene = "Cph",
        dim_red = "UMAP",
        order_cells = TRUE,
        split_conditions = "combined"
    ) +
        ggtitle("Combined Dataset")+
      scale_y_reverse())

    output$vln_expression <- renderPlot(jitter_plot_fun(seurat,
        gene = "Cph",
        celltypes = c(
            "ISC", "EB", "EEP", "dEC", "daEC", "aEC",
            "mEC", "Copper", "LFC", "pEC", "EE"
        ),
        summary_stats = "Include mean",
        split_conditions = "combined"
    ))

    output$annotation_plot_rnai <- renderPlot(annotation_plot_fun(seurat_RNAi,
        dim_red = "UMAP",
        split_conditions = "combined"
    ) +
        ggtitle("Combined Dataset")+scale_x_reverse())

    output$feature_plot_rnai <- renderPlot(feature_plot_fun(seurat_RNAi,
        gene = "esg",
        dim_red = "UMAP",
        order_cells = TRUE,
        split_conditions = "combined"
    ) +
        ggtitle("Combined Dataset")+scale_x_reverse())

    output$vln_expression_rnai <- renderPlot(jitter_plot_fun(seurat_RNAi,
        gene = "esg",
        celltypes = c(
            "ISC", "EB", "EEP", "dEC", "daEC", "aEC",
            "mEC", "Copper", "LFC", "pEC", "EE"
        ),
        summary_stats = "Include mean",
        split_conditions = "combined"
    ))

    ##################
    # NOTCH KO PLOTS #
    ##################
    observeEvent(input$go, {
        dim_red <- input$dim_reduction
        order_cells <- input$order_cells
        celltypes <- input$celltypes_vln_plot
        gene_name <- input$selected_gene
        # if the user selected a flybase gene id, we map it to the symbol
        if (grepl("FBgn\\d+", gene_name)) {
            gene_name <- feature_mapping %>%
                filter(fbid == gene_name) %>%
                pull(symbol)
        }
        summary_stats <- input$summary_stats
        split_condition <- input$split_condition

        if (!split_condition) {
            ###################
            # ANNOTATION PLOT #
            ###################
            output$annotation_plot <- renderPlot(annotation_plot_fun(seurat, dim_red = dim_red, split_conditions = "combined")+scale_y_reverse()) %>%
                bindCache(dim_red, "combined", cache = "session")


            #################
            # Feature PLOT #
            ################
            output$feature_plot <- renderPlot(feature_plot_fun(seurat,
                gene = gene_name,
                dim_red = dim_red,
                order_cells = order_cells,
                split_conditions = "combined"
            )+scale_y_reverse()) %>%
                bindCache(dim_red, gene_name, order_cells, "combined", cache = "session")



            ################
            # Violin PLOT #
            ###############
            output$vln_expression <- renderPlot(jitter_plot_fun(seurat,
                gene = gene_name,
                celltypes = celltypes,
                summary_stats = summary_stats,
                split_conditions = "combined"
            )) %>%
                bindCache(celltypes, gene_name, summary_stats, "combined", cache = "session")
        }

        # if split condition, generate different plots
        if (split_condition) {
            ###################
            # ANNOTATION PLOT #
            ###################
            output$annotation_plot <- renderPlot(annotation_plot_fun(seurat, dim_red = dim_red, split_conditions = "ctrl")+scale_y_reverse()) %>%
                bindCache(dim_red, "ctrl", cache = "session")

            output$annotation_plot_split <- renderPlot(annotation_plot_fun(seurat, dim_red = dim_red, split_conditions = "notch")+scale_y_reverse()) %>%
                bindCache(dim_red, "notch", cache = "session")

            #################
            # Feature PLOT #
            ################
            output$feature_plot <- renderPlot(feature_plot_fun(seurat,
                gene = gene_name,
                dim_red = dim_red,
                order_cells = order_cells,
                split_conditions = "ctrl"
            )+scale_y_reverse()) %>%
                bindCache(dim_red, gene_name, order_cells, "ctrl", cache = "session")

            output$feature_plot_split <- renderPlot(feature_plot_fun(seurat,
                gene = gene_name,
                dim_red = dim_red,
                order_cells = order_cells,
                split_conditions = "notch"
            )+scale_y_reverse()) %>%
                bindCache(dim_red, gene_name, order_cells, "notch", cache = "session")

            ################
            # Violin PLOT #
            ###############
            output$vln_expression <- renderPlot(jitter_plot_fun(seurat,
                gene = gene_name,
                celltypes = celltypes,
                summary_stats = summary_stats,
                split_conditions = "split"
            )) %>%
                bindCache(celltypes, gene_name, summary_stats, "split", cache = "session")
        }
    })

    ##############
    # RNAi PLOTS #
    ##############
    observeEvent(input$go_rnai, {
        dim_red_rnai <- input$dim_reduction_rnai
        order_cells_rnai <- input$order_cells_rnai
        celltypes_rnai <- input$celltypes_vln_plot_rnai
        gene_name_rnai <- input$selected_gene_rnai
        # if the user selected a flybase gene id, we map it to the symbol
        if (grepl("FBgn\\d+", gene_name_rnai)) {
            gene_name_rnai <- feature_mapping %>%
                filter(fbid == gene_name_rnai) %>%
                pull(symbol)
        }
        summary_stats_rnai <- input$summary_stats_rnai
        selected_conditions_rnai <- input$selected_conditions_rnai


        ###################
        # ANNOTATION PLOT #
        ###################
        output$annotation_plot_rnai <- renderPlot(annotation_plot_fun(seurat_RNAi, dim_red = dim_red_rnai, split_conditions = "combined")+scale_x_reverse())


        #################
        # Feature PLOT #
        ################
        output$feature_plot_rnai <- renderPlot(feature_plot_fun(seurat_RNAi,
            gene = gene_name_rnai,
            dim_red = dim_red_rnai,
            order_cells = order_cells_rnai,
            split_conditions = "combined"
        )+scale_x_reverse())



        ################
        # Violin PLOT #
        ###############
        output$vln_expression_rnai <- renderPlot(jitter_plot_fun(seurat_RNAi,
            gene = gene_name_rnai,
            celltypes = celltypes_rnai,
            summary_stats = summary_stats_rnai,
            split_conditions = "combined"
        ))


        if (!is.null(selected_conditions_rnai)) {
            if ("ctrl" %in% selected_conditions_rnai) {
                output$annotation_plot_ctrl_rnai <- renderPlot(annotation_plot_fun(
                    seurat_RNAi,
                    dim_red = dim_red_rnai, split_conditions = "ctrl"
                )+scale_x_reverse())

                output$feature_plot_ctrl_rnai <- renderPlot(feature_plot_fun(seurat_RNAi,
                    gene = gene_name_rnai,
                    dim_red = dim_red_rnai,
                    order_cells = order_cells_rnai,
                    split_conditions = "ctrl"
                )+scale_x_reverse())
            }
            if ("NotchRNAi" %in% selected_conditions_rnai) {
                output$annotation_plot_NotchRNAi_rnai <- renderPlot(annotation_plot_fun(
                    seurat_RNAi,
                    dim_red = dim_red_rnai, split_conditions = "NotchRNAi"
                )+scale_x_reverse())

                output$feature_plot_NotchRNAi_rnai <- renderPlot(feature_plot_fun(seurat_RNAi,
                    gene = gene_name_rnai,
                    dim_red = dim_red_rnai,
                    order_cells = order_cells_rnai,
                    split_conditions = "NotchRNAi"
                )+scale_x_reverse())
            }
            if ("CphUp" %in% selected_conditions_rnai) {
                output$annotation_plot_CphUp_rnai <- renderPlot(annotation_plot_fun(
                    seurat_RNAi,
                    dim_red = dim_red_rnai, split_conditions = "CphUp"
                )+scale_x_reverse())

                output$feature_plot_CphUp_rnai <- renderPlot(feature_plot_fun(seurat_RNAi,
                    gene = gene_name_rnai,
                    dim_red = dim_red_rnai,
                    order_cells = order_cells_rnai,
                    split_conditions = "CphUp"
                )+scale_x_reverse())
            }
            if ("NotchCphRNAi" %in% selected_conditions_rnai) {
                output$annotation_plot_NotchCphRNAi_rnai <- renderPlot(annotation_plot_fun(
                    seurat_RNAi,
                    dim_red = dim_red_rnai, split_conditions = "NotchCphRNAi"
                )+scale_x_reverse())

                output$feature_plot_NotchCphRNAi_rnai <- renderPlot(feature_plot_fun(seurat_RNAi,
                    gene = gene_name_rnai,
                    dim_red = dim_red_rnai,
                    order_cells = order_cells_rnai,
                    split_conditions = "NotchCphRNAi"
                )+scale_x_reverse())
            }
            output$vln_expression_rnai <- renderPlot(jitter_plot_fun_rnai(seurat_RNAi,
                gene = gene_name_rnai,
                celltypes = celltypes_rnai,
                summary_stats = summary_stats_rnai,
                conditions = selected_conditions_rnai
            ))
        }
    })

    #################
    # DownloadCalls #
    #################
    observe({
        if (!input$split_condition) {
            output$download_annotation_plot <- downloadHandler(
                filename = function() {
                    paste0("annotation_", input$dim_reduction, ".pdf")
                },
                content = function(file) {
                    pdf(file, width = 7.31, height = 5.61)
                    print(annotation_plot_fun(seurat,
                        dim_red = input$dim_reduction,
                        split_conditions = "combined"
                    ))
                    dev.off()
                }
            )

            output$download_feature_plot <- downloadHandler(
                filename = function() {
                    paste0("feature_expression_", input$selected_gene, ".pdf")
                },
                content = function(file) {
                    pdf(file, width = 7.31, height = 5.61)
                    print(feature_plot_fun(seurat,
                        gene = input$selected_gene,
                        dim_red = input$dim_reduction,
                        order_cells = input$order_cells,
                        split_conditions = "combined"
                    ))
                    dev.off()
                }
            )

            output$download_vln_expression <- downloadHandler(
                filename = function() {
                    paste0("celltype_expression_", input$selected_gene, ".pdf")
                },
                content = function(file) {
                    p <- jitter_plot_fun(seurat,
                        gene = input$selected_gene,
                        celltypes = input$celltypes_vln_plot,
                        summary_stats = input$summary_stats,
                        split_conditions = "combined"
                    )
                    ggsave(filename = file, width = 7, height = 2.5, plot = p, scale = 1.4)
                }
            )
        } else {
            output$download_annotation_plot <- downloadHandler(
                filename = function() {
                    paste0("annotation_", input$dim_reduction, ".pdf")
                },
                content = function(file) {
                    pdf(file, width = 7.31, height = 5.61)
                    print(annotation_plot_fun(seurat,
                        dim_red = input$dim_reduction,
                        split_conditions = "ctrl"
                    ))
                    dev.off()
                }
            )

            output$download_feature_plot <- downloadHandler(
                filename = function() {
                    paste0("feature_expression_", input$selected_gene, ".pdf")
                },
                content = function(file) {
                    pdf(file, width = 7.31, height = 5.61)
                    print(feature_plot_fun(seurat,
                        gene = input$selected_gene,
                        dim_red = input$dim_reduction,
                        order_cells = input$order_cells,
                        split_conditions = "ctrl"
                    ))
                    dev.off()
                }
            )

            output$download_vln_expression <- downloadHandler(
                filename = function() {
                    paste0("celltype_expression_", input$selected_gene, ".pdf")
                },
                content = function(file) {
                    p <- jitter_plot_fun(seurat,
                        gene = input$selected_gene,
                        celltypes = input$celltypes_vln_plot,
                        summary_stats = input$summary_stats,
                        split_conditions = "split"
                    )
                    ggsave(filename = file, width = 7, height = 2.5, plot = p, scale = 1.4)
                }
            )

            output$download_annotation_plot_split <- downloadHandler(
                filename = function() {
                    paste0("annotation_notch_", input$dim_reduction, ".pdf")
                },
                content = function(file) {
                    pdf(file, width = 7.31, height = 5.61)
                    print(annotation_plot_fun(seurat,
                        dim_red = input$dim_reduction,
                        split_conditions = "notch"
                    ))
                    dev.off()
                }
            )

            output$download_feature_plot_split <- downloadHandler(
                filename = function() {
                    paste0("feature_expression_notch_", input$selected_gene, ".pdf")
                },
                content = function(file) {
                    pdf(file, width = 7.31, height = 5.61)
                    print(feature_plot_fun(seurat,
                        gene = input$selected_gene,
                        dim_red = input$dim_reduction,
                        order_cells = input$order_cells,
                        split_conditions = "notch"
                    ))
                    dev.off()
                }
            )
        }
    })
}

shinyApp(ui, server)
