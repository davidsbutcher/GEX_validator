
library(magrittr)
library(dplyr)
library(ggplot2)
library(shiny)
library(shinyBS)
library(shinyFiles)
library(tippy)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(shinybusy)

# UI ----------------------------------------------------------------------

ui <-
    function(request) {

        fillPage(
            titlePanel("GEX validator"),
            theme = "maglab_theme_old.css",

            useShinyjs(),

            setBackgroundColor(
                color = c("#FFFFFF"),
                gradient = "linear",
                direction = "bottom"
            ),

            tags$head(
                tags$style(HTML("hr {border-top: 1px solid #A9A9A9;}")),
                tags$style(
                    HTML(
                        ".shiny-notification {
                        position:fixed;
                        top: calc(20%);
                        left: calc(40%);
                    }"
                    )
                )
            ),

            add_busy_spinner(
                spin = "fading-circle",
                timeout = 500,
                height = '100px',
                width = '100px'
            ),

            sidebarLayout(

                # Sidebar -----------------------------------------------------------------

                sidebarPanel(
                    width = 2,
                    tabsetPanel(
                        id = "mainpanel",
                        type = "pills",
                        tabPanel(
                            "Upload",
                            br(),
                            shinyDirButton(
                                "GEXresultsdir",
                                "Select GEX results",
                                "Select GEX results folder"
                            ),
                            br(),
                            br(),
                            actionButton(
                                "startValidation",
                                "Start Validation"
                            ),
                            br(),
                            br(),
                            actionButton(
                                "saveSpectra",
                                "Save spectra"
                            ),
                            br(),
                            br(),
                            # actionButton(
                            #     "bookmark",
                            #     "Bookmark state"
                            # ),
                            br(),
                            br(),
                            selectInput(
                                "sequence_i",
                                "SEQ NAME",
                                choices = NULL
                            ),
                            br(),
                            selectInput(
                                "text_size",
                                "TEXT SIZE",
                                choices = c(12,20,28,36)
                            )
                        ),
                        tabPanel(
                            "About",
                            hr(),
                            includeMarkdown("about.md")
                        )
                    )
                ),


                # Main Panel --------------------------------------------------------------

                mainPanel(
                    width = 10,
                    uiOutput("results_test"),
                    textOutput("validation"),
                    h4(
                        splitLayout(
                            splitLayout(
                                cellWidths = 125,
                                "TOTAL SEQS:",
                                textOutput("specnum_total")
                            ),
                            splitLayout(
                                cellWidths = 175,
                                "SEQ NAME:",
                                textOutput("sequence_i")
                            ),
                            splitLayout(
                                cellWidths = c(125,50,25,50),
                                "SPEC NUM:",
                                textOutput("specnum_j"),
                                "/",
                                textOutput("speclength_j")
                            )
                        )
                    ),
                    plotOutput(
                        "plot",
                        width = "1400px",
                        height = "700px"
                    ),
                    actionButton("truepos", "True Positive", icon = icon("check-circle")),
                    actionButton("falsepos", "False Positive", icon = icon("times-circle")),
                    actionButton("goback", "Go Back", icon = icon("backward"))
                ),
                fluid = FALSE
            )
        )

    }


# Server ------------------------------------------------------------------

server <-
    function(input, output, session) {

        # Initial setup -----------------------------------------------------------

        # Establish params to use for shinyFiles input (local only)

        volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())

        shinyDirChoose(
            input = input,
            "GEXresultsdir",
            session = session,
            roots = volumes
        )

        # Disable buttons

        disable("truepos")
        disable("falsepos")
        disable("goback")
        disable("startValidation")
        disable("saveSpectra")


        # Mutable outputs ---------------------------------------------------------

        output$sequence_i <-
            renderText({test_seqs()[iterator$i]})


        output$specnum_j <-
            renderText(
                {
                    validate(
                        need(spec_objects_length(), "")
                    )

                    iterator$j
                }
            )

        output$specnum_total <-
            renderText(
                {
                    validate(
                        need(spec_objects(), "")
                    )

                    length(spec_objects())
                }
            )

        output$speclength_j <-
            renderText({spec_objects_length()})


        # Reactives ---------------------------------------------------------------

        dir_path <-
            reactive(
                {
                    parseDirPath(
                        roots = volumes,
                        input$GEXresultsdir
                    )
                }
            )

        is_valid <-
            reactive(
                {
                    validate(
                        need(dir_path(), "")
                    )

                    files <-
                        fs::dir_ls(
                            dir_path(),
                            type = "file"
                        ) %>%
                        purrr::map_chr(
                            basename
                        )

                    directories <-
                        fs::dir_ls(
                            dir_path(),
                            type = "directory"
                        ) %>%
                        purrr::map_chr(
                            basename
                        )

                    truth_vec <- vector(mode = "logical", length = 3)

                    truth_vec[[1]] <-
                        {
                            all(
                                c("GEX_params.txt", "GEX_timers.csv") %in% files
                            )
                        }

                    truth_vec[[2]] <-
                        {
                            all(
                                c("assigned_fragments_MS2", "IsoDist_MS2", "specZoom_MS2", "XIC_MS2") %in% directories
                            )
                        }

                    truth_vec[[3]] <-
                        {
                            if (truth_vec[[2]] == TRUE) {
                                fs::dir_ls(
                                    fs::path(dir_path(), "specZoom_MS2"),
                                    type = "file"
                                ) %>%
                                    purrr::map_chr(
                                        ~stringr::str_detect(.x, "_specObject_MS2.rds")
                                    ) %>%
                                    as.logical() %>%
                                    any()
                            } else {
                                FALSE
                            }
                        }

                    all(truth_vec)

                }
            )

        test_seqs_path <-
            reactive(
                {
                    validate(
                        need(is_valid(), "")
                    )

                    files <-
                        fs::dir_ls(
                            dir_path(),
                            type = "file"
                        )

                    test_seqs_index <-
                        files %>%
                        purrr::map_chr(
                            ~stringr::str_detect(.x, "target_seqs")
                        ) %>%
                        as.logical() %>%
                        which()

                    files[test_seqs_index]
                }
            )

        test_seqs <-
            reactive(
                {
                    validate(
                        need(is_valid(), "GEX results not valid"),
                        need(spec_objects(), "spec_objects() invalid")
                    )

                    spec_objects() %>%
                        fs::path_file() %>%
                        stringr::str_extract("[^_]+")

                    # Old method for getting test sequence names

                    # readr::read_csv(
                    #     test_seqs_path()
                    # ) %>%
                    #     dplyr::pull(1)
                }
            )

        spec_objects <-
            reactive(
                {
                    validate(
                        need(is_valid(), "")
                    )

                    fs::dir_ls(
                        fs::path(dir_path(), "specZoom_MS2"),
                        type = "file",
                        glob = "*.rds"
                    )
                }
            )

        spec_objects_length <-
            reactive(
                {
                    validate(
                        need(is_valid(), ""),
                        need(spec_objects(), "")
                    )

                    length(readRDS(spec_objects()[[iterator$i]])[[1]])
                }
            )


        reactive_plot <-
            reactive(
                {
                    validate(
                        need(spec_objects(), ""),
                        need(spec_objects_length(), "")
                    )

                    if (spec_objects_length() != 0) {

                        readRDS(spec_objects()[[iterator$i]])[[1]][[iterator$j]]

                    } else {

                        iterator$i <- iterator$i + 1

                        updateSelectInput(
                            session,
                            "sequence_i",
                            choices = test_seqs(),
                            selected = test_seqs()[iterator$i]
                        )

                    }

                }
            )

        plot_theme <-
            reactive(
                {
                    theme(
                        text = element_text(size = input$text_size)
                    )
                }
            )

        # Reactive values ---------------------------------------------------------

        iterator <-
            reactiveValues(
                i = 1,
                j = 1
            )

        validated_plots <-
            reactiveValues()

        invalidated_plots <-
            reactiveValues()

        # Listeners ---------------------------------------------------------------

        listener_upload <-
            reactive(
                {

                    list(
                        input$GEXresultsdir
                    )
                }
            )

        # Observers ---------------------------------------------------------------

        observeEvent(
            input$sequence_i,
            {
                temp_i <- iterator$i
                temp_j <- iterator$j

                iterator$i <- which(test_seqs() == input$sequence_i)
                iterator$j <- 1

                if (spec_objects_length() == 0) {

                    showModal(
                        modalDialog(
                            title = "Invalid sequence",
                            "This sequence has no spectra. Please select a different sequence."
                        )
                    )

                    iterator$i <- temp_i
                    iterator$j <- temp_j

                    updateSelectInput(
                        session,
                        "sequence_i",
                        choices = test_seqs(),
                        selected = test_seqs()[iterator$i]
                    )

                }
            }
        )

        observeEvent(
            listener_upload(),
            {
                # output$results_test <-
                #     renderUI(
                #         {
                #             HTML(
                #                 paste(
                #                     # dir_path(),
                #                     is_valid(),
                #                     # test_seqs_path(),
                #                     # spec_objects(),
                #                     sep = '<br/>'
                #                 )
                #             )
                #         }
                #     )

                iterator$i <- 1
                iterator$j <- 1
                output$results_test <- NULL
                output$plot <- NULL
                validated_plots <- reactiveValues()
                invalidated_plots <- reactiveValues()
                enable("startValidation")

            }
        )

        observeEvent(
            input$startValidation,
            {
                output$validation <- NULL
                output$results_test <- NULL

                disable("startValidation")
                disable("saveSpectra")

                updateSelectInput(
                    session,
                    "sequence_i",
                    choices = test_seqs()
                )

                # Show the first plot

                output$plot <-
                    renderCachedPlot(
                        {
                            reactive_plot()
                        },
                        cacheKeyExpr = paste0(test_seqs()[[iterator$i]], iterator$j, collapse = ""),
                        res = 250
                    )

                # Enable buttons

                enable("truepos")
                enable("falsepos")
                enable("goback")

            }
        )

        observeEvent(
            input$truepos,
            {
                validated_plots[[as.character(iterator$i)]][[as.character(iterator$j)]] <-
                    reactive_plot()

                invalidated_plots[[as.character(iterator$i)]][[as.character(iterator$j)]] <-
                    NULL
            }
        )

        observeEvent(
            input$falsepos,
            {
                validated_plots[[as.character(iterator$i)]][[as.character(iterator$j)]] <-
                    NULL

                invalidated_plots[[as.character(iterator$i)]][[as.character(iterator$j)]] <-
                    reactive_plot()
            }
        )

        observeEvent(
            list(
                input$truepos,
                input$falsepos
            ),
            {
                if (iterator$j < spec_objects_length()) {

                    iterator$j <- iterator$j + 1

                } else if (iterator$i < length(spec_objects())) {

                    iterator$j <- 1
                    iterator$i <- iterator$i + 1

                    updateSelectInput(
                        session,
                        "sequence_i",
                        choices = test_seqs(),
                        selected = test_seqs()[iterator$i]
                    )

                } else {

                    output$results_test <-
                        renderText({"VALIDATION COMPLETE"})

                    disable("truepos")
                    disable("falsepos")
                    disable("goback")
                    enable("saveSpectra")

                }

                output$plot <-
                    renderCachedPlot(
                        {reactive_plot()},
                        cacheKeyExpr = paste0(test_seqs()[[iterator$i]], iterator$j, collapse = ""),
                        res = 250
                    )
            }
        )

        observeEvent(
            input$goback,
            {

                if (iterator$j > 1) {

                    iterator$j <- iterator$j - 1

                } else if (iterator$j == 1 & iterator$i > 1) {

                    iterator$i <- iterator$i - 1
                    iterator$j <- spec_objects_length()

                }
            }
        )

        observeEvent(
            input$saveSpectra,
            {

                output$results_test <-
                    renderText({"Saving results to PDF"})

                # Get raw file name

                GEXparams <-
                    readr::read_lines(
                        file = fs::path(dir_path(), "GEX_params", ext = "txt")
                    )

                fileNameIndex <-
                    GEXparams %>%
                    stringr::str_detect("rawFileName") %>%
                    which()

                rawFile <-
                    GEXparams[fileNameIndex] %>%
                    stringr::str_split_fixed(" = ", n = 2) %>%
                    .[2]

                # Get list of sequence names
                # Might not need this.. causes problems if a sequence
                # has no validated fragments

                # list_names <-
                #     purrr::map_chr(
                #         spec_objects(),
                #         ~basename(.x) %>%
                #             stringr::str_extract("[^_]+")
                #     )

                # Save validated/invalidated plots to lists

                validated_plots_list <-
                    reactiveValuesToList(validated_plots) %>%
                    {.[rev(names(.))]} %>%
                    # purrr::set_names(list_names) %>%
                    purrr::compact()

                invalidated_plots_list <-
                    reactiveValuesToList(invalidated_plots) %>%
                    {.[rev(names(.))]} %>%
                    # purrr::set_names(list_names) %>%
                    purrr::compact()

                # marrange and save plots

                if (length(validated_plots_list) > 0) {

                    valid_spectra_MS2_marrange <-
                        unlist(validated_plots_list, recursive = F) %>%
                        purrr::flatten() %>%
                        gridExtra::marrangeGrob(
                            grobs = .,
                            ncol = 5,
                            nrow = 3,
                            top = rawFile
                        )

                    ggplot2::ggsave(
                        filename =
                            fs::path(
                                dir_path(),
                                paste0(
                                    fs::path_ext_remove(rawFile),
                                    "_validated_fragments",
                                    sep = ""
                                ),
                                ext = "pdf"
                            ),
                        plot = valid_spectra_MS2_marrange,
                        width = 20,
                        height = 12,
                        limitsize = FALSE
                    )

                }

                if (length(invalidated_plots_list) > 0) {

                    invalid_spectra_MS2_marrange <-
                        unlist(invalidated_plots_list, recursive = F) %>%
                        purrr::flatten() %>%
                        gridExtra::marrangeGrob(
                            grobs = .,
                            ncol = 5,
                            nrow = 3,
                            top = rawFile
                        )

                    ggplot2::ggsave(
                        filename =
                            fs::path(
                                dir_path(),
                                paste0(
                                    fs::path_ext_remove(rawFile),
                                    "_invalidated_fragments",
                                    sep = ""
                                ),
                                ext = "pdf"
                            ),
                        plot = invalid_spectra_MS2_marrange,
                        width = 20,
                        height = 12,
                        limitsize = FALSE
                    )

                }

                # Extract data from ggplot labels

                if (length(validated_plots_list) > 0) {

                    label_data_valid <-
                        purrr::map(
                            unlist(validated_plots_list, recursive = F),
                            ~.x[[1]][["layers"]][[2]][["geom_params"]][["label"]]
                        ) %>%
                        stringr::str_split_fixed("\n", n = 4) %>%
                        tibble::as_tibble() %>%
                        dplyr::select(-V4) %>%
                        dplyr::mutate(V3 = as.numeric(stringr::str_extract(V3, "\\d{1,4}"))) %>%
                        purrr::set_names(
                            c("sequence_name", "ion", "charge")
                        )

                }

                if (length(invalidated_plots_list) > 0) {

                    label_data_invalid <-
                        purrr::map(
                            unlist(invalidated_plots_list, recursive = F),
                            ~.x[[1]][["layers"]][[2]][["geom_params"]][["label"]]
                        ) %>%
                        stringr::str_split_fixed("\n", n = 4) %>%
                        tibble::as_tibble() %>%
                        dplyr::select(-V4) %>%
                        dplyr::mutate(V3 = as.numeric(stringr::str_extract(V3, "\\d{1,4}"))) %>%
                        purrr::set_names(
                            c("sequence_name", "ion", "charge")
                        )

                }

                # Load assigned fragments

                ass_frag <-
                    readxl::read_xlsx(
                        fs::path(
                            dir_path(),
                            paste0(fs::path_ext_remove(rawFile), "_assigned_fragments_MS2.xlsx")
                        )
                    )

                # Combine label data with assigned fragments, write spreadsheets

                if (length(validated_plots_list) > 0) {

                    validated_frags_spreadsheet <-
                        dplyr::left_join(
                            label_data_valid,
                            ass_frag
                        ) %>%
                        dplyr::select(
                            raw_filename,
                            dplyr::everything()
                        )

                    writexl::write_xlsx(
                        validated_frags_spreadsheet,
                        fs::path(
                            dir_path(),
                            paste0(fs::path_ext_remove(rawFile), "_validated_fragments_MS2.xlsx")
                        )
                    )

                }

                if (length(invalidated_plots_list) > 0) {

                    invalidated_frags_spreadsheet <-
                        dplyr::left_join(
                            label_data_invalid,
                            ass_frag
                        ) %>%
                        dplyr::select(
                            raw_filename,
                            dplyr::everything()
                        )

                    writexl::write_xlsx(
                        invalidated_frags_spreadsheet,
                        fs::path(
                            dir_path(),
                            paste0(fs::path_ext_remove(rawFile), "_invalidated_fragments_MS2.xlsx")
                        )
                    )

                }

            }
        )

        observeEvent(
            input$bookmark,
            {session$doBookmark()}
        )


        # Bookmarking -------------------------------------------------------------

        # Save extra values in state$values when we bookmark
        onBookmark(
            function(state) {
                state$values$iterator_i <- iterator$i
                state$values$iterator_j <- iterator$j
                state$values$validated_plots <- validated_plots
                state$values$invalidated_plots <- invalidated_plots
            }
        )

        # Read values from state$values when we restore
        onRestore(
            function(state) {
                iterator$i <- state$values$iterator_i
                iterator$j <- state$values$iterator_j
                validated_plots <- state$values$validated_plots
                invalidated_plots <- state$values$invalidated_plots
            }
        )


    }


# Run the application

shinyApp(ui = ui, server = server, enableBookmarking = "server")
