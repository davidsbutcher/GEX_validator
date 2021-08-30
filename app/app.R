
library(magrittr)
library(dplyr)
library(ggplot2)
library(DT)
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
                timeout = 250,
                height = '75px',
                width = '75px'
            ),

            sidebarLayout(

                # Sidebar -----------------------------------------------------------------

                sidebarPanel(
                    width = 2,
                    tabsetPanel(
                        id = "sidebarpanel",
                        type = "pills",
                        tabPanel(
                            "Load data",
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
                                "Start validation"
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
                            br()
                            # selectInput(
                            #     "text_size",
                            #     "TEXT SIZE",
                            #     choices = c(12,20,28,36)
                            # )
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
                    tabsetPanel(
                        id = "mainpanel",
                        type = "tabs",
                        tabPanel(
                            "Plot",
                            uiOutput("results_test"),
                            h4(
                                splitLayout(
                                    splitLayout(
                                        textOutput("validation")
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
                            h5(
                                style = "color:red;text-align:left;",
                                splitLayout(
                                    splitLayout(
                                        cellWidths = 175,
                                        "NEXT SEQ:",
                                        textOutput("nextseq_seqname")
                                    )
                                ),
                                splitLayout(
                                    splitLayout(
                                        cellWidths = 175,
                                        "ION:",
                                        textOutput("nextseq_ion")
                                    )
                                ),
                                splitLayout(
                                    splitLayout(
                                        cellWidths = 175,
                                        "CHARGE:",
                                        textOutput("nextseq_charge")
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
                            # br(),
                            # br(),
                            # actionButton("validateRemaining", "VALIDATE REMAINING", icon = icon("skull")),
                            # actionButton("invalidateRemaining", "INVALIDATE REMAINING", icon = icon("skull"))
                        ),
                        tabPanel(
                            "Validated",
                            div(
                                DTOutput("validtable", height = "700px"),
                                style = "width: 1200px; height: 700px; overflow-y: scroll;"
                            )
                        )
                    )
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

        output$nextseq_seqname <-
            renderText({reactive_next_plot_info()$sequence_name})

        output$nextseq_ion <-
            renderText({reactive_next_plot_info()$ion})

        output$nextseq_charge <-
            renderText({reactive_next_plot_info()$charge})

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

        output$validation <-
            renderText(
                {
                    validate(
                        need(is_valid(), "GEX results not valid"),
                        need(spec_objects(), "spec_objects() invalid")
                    )


                    if (length(reactiveValuesToList(complete_seqs)) == length(spec_objects())) {

                        # disable("truepos")
                        # disable("falsepos")
                        # disable("goback")
                        # enable("saveSpectra")
                        "ALL SEQUENCES COMPLETE"

                    } else {

                        paste0(length(reactiveValuesToList(complete_seqs)), "/", length(spec_objects()), " SEQUENCES COMPLETE", collapse = "")

                    }
                }
            )

        output$validtable <-
            renderDT(
                {
                    validate(
                        need(test_seqs(), "")
                    )

                    val_comp <- vector(mode = "logical", length = length(test_seqs()))
                    comp_seqs <- unlist(reactiveValuesToList(complete_seqs))

                    val_comp[which(comp_seqs == TRUE)] <- TRUE

                    tibble::tibble(
                        `Sequence Name` = test_seqs(),
                        `Validation Complete` = val_comp
                    )
                },
                options =
                    list(
                        pageLength = 15,
                        lengthChange = FALSE
                    )
            )

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
                                c("GEX_params.txt") %in% files
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

                        complete_seqs[[as.character(iterator$i)]] <- TRUE

                        updateSelectInput(
                            session,
                            "sequence_i",
                            choices = test_seqs(),
                            selected = test_seqs()[iterator$i]
                        )

                    }

                }
            )

        reactive_next_plot_info <-
            reactive(
                {
                    validate(
                        need(spec_objects(), ""),
                        need(spec_objects_length(), "")
                    )

                    if (spec_objects_length() != 0 && spec_objects_length() > iterator$j) {

                        readRDS(spec_objects()[[iterator$i]])[[1]][[iterator$j+1]] %>%
                            {.[[1]][["layers"]][[2]][["computed_geom_params"]][["label"]]} %>%
                            stringr::str_split_fixed("\n", n = 5) %>%
                            tibble::as_tibble(.name_repair = "unique") %>%
                            dplyr::select(1, 2, 4) %>%
                            purrr::set_names(
                                c("sequence_name", "ion", "charge")
                            ) %>%
                            dplyr::mutate(charge = as.numeric(stringr::str_extract(charge, "\\d{1,4}")))

                    } else {

                        tibble::tibble(
                            "sequence_name" = NA,
                            "ion" = NA,
                            "charge" = NA
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

        checked_plots_number <-
            reactive(
                {
                    valplots <-
                        validated_plots[[as.character(iterator$i)]] %>%
                        purrr::map(
                            ~!is.null(.x)
                        ) %>%
                        unlist() %>%
                        sum()

                    invalplots <-
                        invalidated_plots[[as.character(iterator$i)]] %>%
                        purrr::map(
                            ~!is.null(.x)
                        ) %>%
                        unlist() %>%
                        sum()

                    valplots + invalplots

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

        complete_seqs <-
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
            input$bookmark,
            {session$doBookmark()}
        )

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
                complete_seqs <- reactiveValues()
                enable("startValidation")

            }
        )

        observeEvent(
            input$startValidation,
            {
                # output$validation <- NULL
                output$results_test <- NULL

                disable("startValidation")
                enable("saveSpectra")

                updateSelectInput(
                    session,
                    "sequence_i",
                    choices = test_seqs()
                )

                if (spec_objects_length() == checked_plots_number()) {
                    complete_seqs[[as.character(iterator$i)]] <- TRUE
                }

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
                # Remove the output$results_test value if truepos or falsepos
                # are pressed

                output$results_test <- NULL

                if (spec_objects_length() == checked_plots_number()) {
                    complete_seqs[[as.character(iterator$i)]] <- TRUE
                }

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

                    # output$results_test <-
                    #     renderText({"VALIDATION COMPLETE"})
                    #
                    # disable("truepos")
                    # disable("falsepos")
                    # disable("goback")
                    # enable("saveSpectra")

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
                # Remove the output$results_test value if goback
                # is pressed

                output$results_test <- NULL

                if (iterator$j > 1) {

                    iterator$j <- iterator$j - 1

                } else if (iterator$j == 1 & iterator$i > 1) {

                    iterator$i <- iterator$i - 1
                    iterator$j <- spec_objects_length()

                }
            }
        )

        observeEvent(
            input$validateRemaining,
            {
                ask_confirmation(
                    "confirmValRemain",
                    title = "Validate remaining fragments",
                    text = "Are you sure you wish to validate all remaining fragments?",
                    type = "question",
                    btn_labels = c("Cancel", "Confirm"),
                    btn_colors = NULL,
                    closeOnClickOutside = FALSE,
                    showCloseButton = FALSE,
                    html = FALSE,
                    session = shiny::getDefaultReactiveDomain()
                )
            }
        )

        observeEvent(
            input$confirmValRemain,
            {
                # This code is nonfunctional. Shelved until I can figure it out. - DSB

                length_spec_objects_temp <- length(spec_objects())
                spec_objects_length_temp <- spec_objects_length()
                spec_objects_temp <- spec_objects()

                for (i in seq_along(length_spec_objects_temp)) {

                    for (j in seq_along(spec_objects_length_temp)) {

                        validated_plots[[as.character(i)]][[as.character(j)]] <-
                            readRDS(spec_objects_temp[[i]])[[1]][[j]]

                    }

                    complete_seqs[[as.character(i)]] <- TRUE

                }
            }
        )

        observeEvent(
            input$invalidateRemaining,
            {
            }
        )

        observeEvent(
            input$saveSpectra,
            {

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

                # Load assigned fragments

                ass_frag <-
                    readxl::read_xlsx(
                        fs::path(
                            dir_path(),
                            paste0(fs::path_ext_remove(rawFile), "_assigned_fragments_MS2.xlsx")
                        )
                    )

                # Save validated plots to list

                validated_plots_list <-
                    reactiveValuesToList(validated_plots) %>%
                    {.[rev(names(.))]} %>%
                    # purrr::set_names(list_names) %>%
                    purrr::compact()


                if (length(validated_plots_list) > 0) {

                    # Extract data from ggplot labels

                    label_data_valid <-
                        purrr::map(
                            unlist(validated_plots_list, recursive = F),
                            ~.x[[1]][["layers"]][[2]][["computed_geom_params"]][["label"]]
                        ) %>%
                        stringr::str_split_fixed("\n", n = 5) %>%
                        tibble::as_tibble() %>%
                        dplyr::select(V1, V2, V4) %>%
                        dplyr::mutate(V4 = as.numeric(stringr::str_extract(V4, "\\d{1,4}"))) %>%
                        purrr::set_names(
                            c("sequence_name", "ion", "charge")
                        )

                    # Combine label data with assigned fragments, write spreadsheets

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

                    # marrange and save validated plots

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

                    # Remove validated plots to save space

                    rm(validated_plots_list)
                    rm(valid_spectra_MS2_marrange)

                }

                # Save invalidated plots to list

                invalidated_plots_list <-
                    reactiveValuesToList(invalidated_plots) %>%
                    {.[rev(names(.))]} %>%
                    # purrr::set_names(list_names) %>%
                    purrr::compact()

                # marrange and save invalidated plots

                if (length(invalidated_plots_list) > 0) {

                    # Extract data from ggplot labels

                    label_data_invalid <-
                        purrr::map(
                            unlist(invalidated_plots_list, recursive = F),
                            ~.x[[1]][["layers"]][[2]][["computed_geom_params"]][["label"]]
                        ) %>%
                        stringr::str_split_fixed("\n", n = 5) %>%
                        tibble::as_tibble() %>%
                        dplyr::select(V1, V2, V4) %>%
                        dplyr::mutate(V4 = as.numeric(stringr::str_extract(V4, "\\d{1,4}"))) %>%
                        purrr::set_names(
                            c("sequence_name", "ion", "charge")
                        )

                    # Combine label data with assigned fragments, write spreadsheets

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

                    # marrange and save invalidated plots

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

                    # Remove validated plots to save space

                    rm(invalidated_plots_list)
                    rm(invalid_spectra_MS2_marrange)

                }

                # Let the user know it worked

                output$results_test <-
                    renderText({"RESULTS SAVED"})

            }
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
