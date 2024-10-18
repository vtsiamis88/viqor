#########################
# VIQoR - Visually Supervised Protein Inference and Protein Quantification
#########################
#Developer and maintainer: Vasileios Tsiamis (vasileios@bmb.sdu.dk)
#########################

#Dependencies
packages <- c("shiny", "shinydashboard", "DT", "protr", "shinyjs", "igraph", "networkD3", 
              "plotly", "shinyBS", "shinycssloaders", "ggplot2", "reshape2", 
              "webshot2", "htmlwidgets", "V8", "seqinr", "dplyr", "seqinr", "heatmaply")

bioconductor <- c("Biostrings", "preprocessCore")

#Install packages
to.install <- setdiff(packages, rownames(installed.packages()))
to.install.b <- setdiff(bioconductor, rownames(installed.packages()))

if(length(to.install) > 0){
  
  install.packages(setdiff(to.install, "webshot2"))
  
}

if(length(which(to.install == "webshot2"))){
  
  install.packages("remotes")
  remotes::install_github("rstudio/webshot2")
  
}

if(length(to.install.b) > 0){
  
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    
    install.packages("BiocManager")
    
  }
    
  BiocManager::install(to.install.b)
  
}

library(shiny)
library(shinydashboard)
library(DT)
library(protr)
library(shinyjs)
library(igraph)
library(networkD3)
library(plotly)
library(shinyBS)
library(shinycssloaders)
library(parallel)
library(ggplot2)
library(reshape2)
library(webshot)
library(htmlwidgets)
library(V8) #shinyjs dependency
library(seqinr)

source("Functions.R")

#USER INTERFACE SIDE
ui <- shinydashboard::dashboardPage(
  
  title = "VIQoR",
  skin = "black",

  shinydashboard::dashboardHeader(title = shiny::tags$div(shiny::tags$img(src = "logo.svg", width = 220)), 
                                  titleWidth = 300,
                                  shiny::tags$li(class = "dropdown", shiny::tags$style(shiny::HTML(".text-info {color:#2A7F69;}"))),
                                  shinydashboard::dropdownMenuOutput("notification_dropdown_menu")), #dasboardHeader end

  shinydashboard::dashboardSidebar( width = 300, useShinyjs(), uiOutput("ui_sidebar")), #dashboardSidebar end

  shinydashboard::dashboardBody(shinyjs::useShinyjs(),
                                shiny::tags$head(shiny::tags$link(rel = "icon", type = "image/png", sizes = "12x12", href = "v.png")),
                                shinyjs::extendShinyjs(text = "shinyjs.collapse = function(boxid) {$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();}", functions=c("collapse")),
                                shiny::uiOutput("ui_body")) #dashboardBody end

)

#SERVER SIDE
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=170*1024^2)
  
  ########################### APPLICATION's GUI ########################### 

  shinydashboard::updateTabItems(session, inputId = "sidebar_menu", selected = "input")
  
  ########################### START dashboardSidebar UI ###########################
  output$ui_sidebar <- shiny::renderUI ({
    
    shinydashboard::sidebarMenu(id = "sidebar_menu",
      
      #Tab 1: input
      shinydashboard::menuItem("Input", tabName = "input", icon = icon("far fa-file"), selected = TRUE),
        shiny::conditionalPanel(condition = "input.sidebar_menu == 'input'",
                                shiny::tags$head(shiny::tags$style(type = "text/css", ".shiny-input-container {padding-top: 0px !important;}")), #remove margin top from shiny objects
                                shiny::tags$hr(id = "input_line1", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em"),
                                shiny::tags$hr(id = "input_line2", style = "border-color: #918785; margin-top: 1px; margin-bottom: 0em"),
                                shiny::tags$h4(id = "input_text1", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("CSV FILE")),
                                shiny::tags$hr(id = "input_line3", style = "border-color: #918785; margin-top: 0em; margin-bottom: 1px"),
                                shiny::tags$hr(id = "input_line4", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em"),
                                shinyjs::hidden(shiny::tags$div(id = "input_cui_csv", shiny::uiOutput("input_file_name_csv"))),
                                shiny::tags$h5(id = "input_empty1", ""),
                                shiny::tags$p(shiny::tags$h5(id = "input_text2", style="display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Choose CSV file:")),
                                              shiny::tags$h6(id = "input_q_icon1", style = "display:inline;", icon("question-circle"))),
                                shinyBS::bsTooltip(id = "input_q_icon1", title = "Imported CSV file must contain the peptide sequences in the first column. At least 4 quantitative columns needed in the imported file. Quantitative columns must be at the end of the table and their identifiers should be unique. Columns between peptide sequences and quantitative measurements will be excluded.", placement = "right", trigger = "hover", options = list(container = "body")),
                                shiny::tags$style(shiny::HTML(".progress-bar {background-color: #2A9C7F}")),
                                shiny::fileInput(inputId = "input_csv", label = NULL, accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                                shiny::tags$div(style="display: inline-block;", shiny::actionLink(inputId = "input_example_csv_link", label = "Example dataset")),
                                                shiny::tags$h6(id = "input_q_icon1_5", style = "display: inline-block;", icon("question-circle"),
                                                shinyBS::bsTooltip(id = "input_q_icon1_5", title = "Click the action link to load an example dataset. The example CSV file is a peptide report from a label-free MS study of hybrid proteome samples with known protein concentrations of 1:1, 10:1 and 1:10 for human, yeast and E.coli respectively. The dataset is taken from: A multicenter study benchmarks software tools for label-free proteome quantification, 2016." , placement = "right", trigger = "hover", options = list(container = "body"))),
                                shiny::tags$h5(id = "input_text3", shiny::HTML("&nbsp;&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Choose CSV file settings:")),
                                shiny::tags$p(shiny::tags$h5(id = "input_text4", style = "display:inline;", shiny::HTML('&nbsp;&nbsp;&nbsp;&nbsp;'), shiny::tags$b("Header:")),
                                              shiny::tags$h6(id = "input_q_icon2", style = "display:inline;", icon("question-circle"))),
                                shinyBS::bsTooltip(id = "input_q_icon2", title = "Indicates whether the CSV file contains the column identifiers as its first row." , placement = "right", trigger = "hover", options = list(container = "body")),
                                shiny::tags$style("input[type='radio']:checked {filter: hue-rotate(320deg);}"),
                                shiny::tags$style("input[type='checkbox']:checked {filter: hue-rotate(320deg);}"),
                                shiny::checkboxInput(inputId = 'input_header_csv', label = "Header?", value = TRUE),
                                shiny::tags$p(shiny::tags$h5(id = "input_text5", style = "display:inline;", shiny::HTML('&nbsp;&nbsp;&nbsp;&nbsp;'), shiny::tags$b("Separator:")),
                                              shiny::tags$h6(id = "input_q_icon3", style = "display:inline;", icon("question-circle"))),
                                shinyBS::bsTooltip(id = "input_q_icon3", title = "The character that separates the values of each row in the CSV file.", placement = "right", trigger = "hover", options = list(container = "body")),
                                shiny::radioButtons(inputId = "input_sep_csv", label = NULL, choices = c(";", ".", ",", ":", "tab"), selected = ",", inline = TRUE),
                                shiny::tags$p(shiny::tags$h5(id = "input_text6", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Decimal separator:")),
                                              shiny::tags$h6(id = "input_q_icon4", style = "display:inline;", icon("question-circle"))),
                                shinyBS::bsTooltip("input_q_icon4", title = "The character that is used in the CSV file as decimal separator.", placement = "right", trigger = "hover", options = list(container = "body")),
                                shiny::radioButtons(inputId = "input_dec_csv", label = NULL, choices = c(".", ","), inline = TRUE),
                                shiny::tags$p(shiny::tags$h5(id = "input_text7", style = "display:inline;", shiny::HTML('&nbsp;&nbsp;&nbsp;&nbsp;'), shiny::tags$b("Data tranformation:")),
                                              shiny::tags$h6(id = "input_q_icon5", style = "display:inline", icon("question-circle"))),
                                shinyBS::bsTooltip(id = "input_q_icon5", title = "Indicates whether the peptide quantitative data in CSV file are log2-transformed.", placement = "right", trigger = "hover", options = list(container = "body")),
                                shiny::checkboxInput(inputId = "input_log_transformed", label = "Data log transformed?", value = FALSE),
                                shiny::tags$hr(id = "input_line5", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em"),
                                shiny::tags$hr(id = "input_line6", style = "border-color: #918785; margin-top: 1px; margin-bottom: 0em"),
                                shiny::tags$h4(id = "input_text8", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("FASTA FILE")),
                                shiny::tags$hr(id = "input_line7", style = "border-color: #918785; margin-top: 0em; margin-bottom: 1px"),
                                shiny::tags$hr(id = "input_line8", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em"),
                                shinyjs::hidden(shiny::tags$div(id = "input_cui_fasta", shiny::uiOutput("input_file_name_fasta"))),
                                shiny::tags$p(shiny::tags$h5(id = "input_text9", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Choose FASTA file:")),
                                              shiny::tags$h6(id = "input_q_icon6", style = "display:inline;", icon("question-circle"))),
                                shinyBS::bsTooltip(id = "input_q_icon6", title = "The FASTA file should contain the sequences of whole proteomes or a mixture or proteomes or a list of specific proteins.", placement = "right", trigger = "hover", options = list(container = "body")),
                                shiny::fileInput(inputId = "input_fasta", label = NULL, accept = c("text/csv", "text/comma-separated-values,text/plain", ".fasta")),
                                shiny::uiOutput("input_action")),
      #Tab 2: conditions_n_samples
      shinydashboard::menuItem("Conditions and Samples", tabName = "conditions_n_samples", icon = icon("cog", lib = "glyphicon")),
        shiny::conditionalPanel(condition = "input.sidebar_menu == 'conditions_n_samples'",
                                shinyjs::hidden(shiny::tags$div(id = "cns_header",
                                  shiny::tags$hr(id = "cns_line1", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em"),
                                  shiny::tags$hr(id = "cns_line2", style = "border-color: #918785; margin-top: 1px; margin-bottom: 0em"),
                                  shiny::tags$h4(id = "cns_text1", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("CONDITIONS AND SAMPLES")),
                                  shiny::tags$hr(id = "cns_line3", style = "border-color: #918785; margin-top: 0em; margin-bottom: 1px"),
                                  shiny::tags$hr(id = "cns_line4", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em")
                                )),
                                shinyjs::hidden(shiny::tags$div(id = "cns_parameters", shiny::uiOutput("cns_parameters_cui"))),
                                shinyjs::hidden(shiny::tags$div(id = "cns_blank1", shiny::tags$h5(id = "cns_empty1", ""))),
                                shiny::tags$style(shiny::HTML(".irs--shiny .irs-max, .irs--shiny .irs-min {color: #cccccc; font-size: 9px;}")),
                                shiny::tags$style(shiny::HTML(".irs--shiny .irs-single {background-color: #2A9C7F; font-size: 10px;}")),
                                shiny::tags$style(shiny::HTML(".irs--shiny .irs-bar {border-top: 1px solid #2A7F69; border-bottom: 1px solid #2A7F69; background: #2A9C7F;}")),
                                shinyjs::hidden(shiny::tags$div(id = "cns_tip1", shiny::tags$p(shiny::tags$h5(id = "cns_text2", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Number of quantitative columns:")),
                                                                                               shiny::tags$h6(id = "cns_q_icon1", style = "display:inline;", icon("question-circle"))),
                                                                shinyBS::bsTooltip(id = "cns_q_icon1", title = "Select the total number of columns in the CSV file that contain peptide quantitative data.The minimum number is set to 4, since it is a requirement for the fast-FARMS summarization method.", placement = "right", trigger = "hover", options = list(container = "body")))),
                                shinyjs::hidden(shiny::tags$div(id = "cns_quantitative_columns_slider", shiny::sliderInput(inputId = "cns_total_quantitative_columns", label = NULL, min = 4, max = 4, step = 1, value = 4, ticks = F))),
                                shinyjs::hidden(shiny::tags$div(id = "cns_tip2", shiny::tags$p(shiny::tags$h5(id = "cns_text3", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Number of conditions:")),
                                                                                               shiny::tags$h6(id = "cns_q_icon2", style = "display:inline;", icon("question-circle"))),
                                                                shinyBS::bsTooltip(id = "cns_q_icon2", title = "Select the total number of conditions. All conditions should have exactly the same number of replicates.", placement = "right", trigger = "hover", options = list(container = "body")))),
                                shinyjs::hidden(shiny::tags$div(id = "cns_conditions_slider", shiny::sliderInput(inputId = "cns_total_conditions", label = NULL, min = 1, max = 1, step = 1, value = 1, ticks = F))),
                                shinyjs::hidden(shiny::tags$div(id = "cns_tip3", shiny::tags$p(shiny::tags$h5(id = "cns_text4", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Reference condition:")),
                                                                                               shiny::tags$h6(id = "cns_q_icon3", style = "display:inline;", icon("question-circle"))),
                                                                shinyBS::bsTooltip(id = "cns_q_icon3", title = "Check the box if there is a reference condition and specify that condition in the appearing slider.", placement = "right", trigger = "hover", options = list(container = "body")))),
                                shinyjs::hidden(shiny::tags$div(id = "cns_reference_checkbox", shiny::checkboxInput(inputId = "cns_q_reference", label = "Reference condition?", value = FALSE))),
                                shiny::uiOutput("cns_reference_ui"),
                                shinyjs::hidden(shiny::tags$div(id = "cns_tip5", shiny::tags$p(shiny::tags$h5(id = "cns_text6", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Column arrangement:")),
                                                                                               shiny::tags$h6(id = "cns_q_icon5", style = "display:inline;", icon("question-circle"))),
                                                                shinyBS::bsTooltip(id = "cns_q_icon5", title = "Check the box if the quantitative columns of the imported CSV file follow an arrangment based on run indices. For example: C1_R1, C2_R1, C3_1, C1_R2, C2_R2, C3_2, C1_R3, C2_R3, C3_3.", placement = "right", trigger = "hover", options = list(container = "body")))),
                                shinyjs::hidden(shiny::tags$div(id = "cns_arrangement_checkbox", shiny::checkboxInput(inputId = "cns_arrangement", label = "Run based arrangement?", value = FALSE))),
                                shinyjs::hidden(shiny::tags$div(id = "cns_tip4", shiny::tags$p(shiny::tags$h5(id = "cns_text5", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Maximum NA values per peptide:")),
                                                                                               shiny::tags$h6(id = "cns_q_icon4", style = "display:inline;", icon("question-circle"))),
                                                                shinyBS::bsTooltip(id = "cns_q_icon4", title = "Slider to choose the maximum amount of missing values allowed per peptide. Peptides with more missing values than the selected number are removed.", placement = "right", trigger = "hover", options = list(container = "body")))),
                                shinyjs::hidden(shiny::tags$div(id = "cns_na_slider", shiny::sliderInput(inputId = "cns_maximum_na", label = NULL, min = 1, max = 1, step = 1, value = 1, ticks = F))),
                                shiny::uiOutput("cns_action")),
      #Tab 3: modifications_n_normalization
      shinydashboard::menuItem("Modifications and Normalization", tabName = "modifications_n_normalization", icon = icon("cog", lib = "glyphicon")),
        shiny::conditionalPanel(condition = "input.sidebar_menu == 'modifications_n_normalization'",
                                shinyjs::hidden(shiny::tags$div(id = "mnn_header1",
                                  shiny::tags$hr(id = "mnn_line1", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em"),
                                  shiny::tags$hr(id = "mnn_line2", style = "border-color: #918785; margin-top: 1px; margin-bottom: 0em"),
                                  shiny::tags$h4(id = "mnn_text1", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("MODIFICATIONS")),
                                  shiny::tags$hr(id = "mnn_line3", style = "border-color: #918785; margin-top: 0em; margin-bottom: 1px"),
                                  shiny::tags$hr(id = "mnn_line4", style="border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em")
                                )),
                                shinyjs::hidden(shiny::tags$div(id = "mnn_blank1", shiny::tags$h5(id = "mnn_empty1", ""))),
                                shiny::uiOutput("mnn_modification_types_tip"),
                                shinyjs::hidden(shiny::tags$div(id = "mnn_modification_types_checkbox", shiny::checkboxGroupInput(inputId = "mnn_modification_types", label = NULL, selected = NULL))),
                                shinyjs::hidden(shiny::tags$div(id = "mnn_blank2", shiny::tags$h5(id = "mnn_empty2"))),
                                shiny::uiOutput("mnn_modifications_action"),
                                shiny::uiOutput("mnn_parameters_cui1"),
                                shinyjs::hidden(shiny::tags$div(id = "mnn_header2",
                                  shiny::tags$hr(id = "mnn_line5", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em"),
                                  shiny::tags$hr(id = "mnn_line6", style = "border-color: #918785; margin-top: 1px; margin-bottom: 0em"),
                                  shiny::tags$h4(id = "mnn_text3", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("NORMALIZATION")),
                                  shiny::tags$hr(id = "mnn_line7", style = "border-color: #918785; margin-top: 0em; margin-bottom: 1px"),
                                  shiny::tags$hr(id = "mnn_line8", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em")
                                )),
                                shinyjs::hidden(shiny::tags$div(id = "mnn_tip2", shiny::tags$p(shiny::tags$h5(id = "mnn_text4", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Normalization method:")),
                                                                                               shiny::tags$h6(id = "mnn_q_icon2", style = "display:inline;", icon("question-circle"))),
                                                                shinyBS::bsTooltip(id = "mnn_q_icon2", title = "Choose the log2 zero center normalization method that should be applied on the peptide abundances. In case of modifications, only unmodified peptide abundances will be normalized.", placement = "right", trigger = "hover", options = list(container = "body")))),
                                shinyjs::hidden(shiny::tags$div(id = "mnn_normalization_radiobuttons", shiny::radioButtons(inputId = "mnn_normalization_method", label = NULL, choices = c("average", "median", "quantile"), selected = "median", inline = TRUE))),
                                shiny::uiOutput("mnn_parameters_cui2"),
                                shiny::uiOutput("mnn_action")),
      #Tab 4: protein_inference
      shinydashboard::menuItem("Protein Inference", tabName = "protein_inference", icon = icon("indent-left", lib = "glyphicon")),
        shiny::conditionalPanel(condition = "input.sidebar_menu == 'protein_inference'",
                                shinyjs::hidden(shiny::tags$div(id = "pi_header1",
                                  shiny::tags$hr(id = "pi_line1", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em"),
                                  shiny::tags$hr(id = "pi_line2", style = "border-color: #918785; margin-top: 1px; margin-bottom: 0em"),
                                  shiny::tags$h4(id = "pi_text1", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("PROTEIN INFERENCE")),
                                  shiny::tags$hr(id = "pi_line3", style = "border-color: #918785; margin-top: 0em; margin-bottom: 1px"),
                                  shiny::tags$hr(id = "pi_line4", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em")
                                )),
                                shinyjs::hidden(shiny::tags$div(id = "pi_blank1", shiny::tags$h5(id = "pi_empty1", ""))),
                                shiny::uiOutput("pi_parameters_cui"),
                                shinyjs::hidden(shiny::tags$div(id = "pi_tip1", shiny::tags$p(shiny::tags$h5(id = "pi_text2", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Inference method:")),
                                                                                              shiny::tags$h6(id = "pi_q_icon1", style = "display:inline;", icon("question-circle"))),
                                                                shinyBS::bsTooltip(id = "pi_q_icon1", title = "Select the protein inference method. Soft and Strict methods use the principle of parsimony. The difference between soft and strict method is the way they handle degenerate peptides. No option directly assigns peptides to the proteins they identify. Strict method creates minimal protein group lists that explains the presense of the peptides uniquely. Soft method creates minimal protein group lists that explains the presense of the peptides at least once. No method creates maximal protein group lists that explain the presense of the peptides at least once.", placement = "right", trigger = "hover", options = list(container = "body")))),
                                shinyjs::hidden(shiny::tags$div(id = "pi_parsimony_radiobuttons", shiny::radioButtons(inputId = "pi_parsimony_type", label = NULL, choices = c("soft", "strict", "no"), selected = character(0), inline = TRUE))),
                                shiny::uiOutput("pi_parsimony_action"),
                                shiny::uiOutput("pi_action")),
      #Tab 5: protein_quantification
      shinydashboard::menuItem("Protein Quantification", tabName = "protein_quantification", icon = icon("bar-chart-o")),
        shiny::conditionalPanel(condition = "input.sidebar_menu == 'protein_quantification'",
                                shinyjs::hidden(shiny::tags$div(id = "pq_header1",
                                  shiny::tags$hr(id = "pq_line1", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em"),
                                  shiny::tags$hr(id = "pq_line2", style = "border-color: #918785; margin-top: 1px; margin-bottom: 0em"),
                                  shiny::tags$h4(id = "pq_text1", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("PROTEIN QUANTIFICATION")),
                                  shiny::tags$hr(id = "pq_line3", style = "border-color: #918785; margin-top: 0em; margin-bottom: 1px"),
                                  shiny::tags$hr(id = "pq_line4", style = "border-color: #D4CCCA; margin-top: 0em; margin-bottom: 0em")
                                )),
                                shinyjs::hidden(shiny::tags$div(id = "pq_blank1", shiny::tags$h5(id = "pq_empty1", ""))),
                                shinyjs::hidden(shiny::tags$div(id = "pq_tip1", shiny::tags$p(shiny::tags$h5(id = "pq_text2", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Fast-FARMS settings:")),
                                                                                              shiny::tags$h6(id = "pq_q_icon1", style = "display:inline;", icon("question-circle"))),
                                                                shinyBS::bsTooltip(id = "pq_q_icon1", title = "Weight hyperparameter value in the range of [0,1] and determines the influence of the prior distribution. Mu hyperparameter value allows to quantify different aspects of potential prior knowledge. Values near zero assumes that most peptides do not contain a signal.", placement = "right", trigger = "hover", options = list(container = "body")))),
                                shinyjs::hidden(shiny::tags$div(id ="pq_farms_weight_slider", shiny::sliderInput(inputId = "pq_farms_weight", label = "FARMS weight parameter:", min = 0, max = 1, step = 0.1, value = 0.1, ticks = F))),
                                shinyjs::hidden(shiny::tags$div(id = "pq_farms_mu_slider", shiny::sliderInput(inputId = "pq_farms_mu", label = "FARMS mu parameter:", min = 0, max = 1, step = 0.1, value = 0.1, ticks = F))),
                                shiny::uiOutput("pq_probes_rescale_cui"),
                                shiny::uiOutput("pq_quantification_action"),
                                shiny::uiOutput("pq_parameters_cui")
                              )
    )
  })
  ########################### END dashboardSidebar UI ###########################

  ########################### START dashboardBody UI ###########################
  output$ui_body <- shiny::renderUI({

    shinydashboard::tabItems(
      #Tab 1: input
      shinydashboard::tabItem(tabName = "input",
                              shiny::tags$style(shiny::HTML(".box-header h3 {font-weight: bold;}")),
                              shiny::tags$style(shiny::HTML(".box.box-solid.box-primary>.box-header {color:#fff;background:#2A9C7F} .box.box-solid.box-primary{border-bottom-color:#2A7F69;border-left-color:#2A7F69;border-right-color:#2A7F69;border-top-color:#2A7F69;}")),
                              shiny::fluidRow(shinydashboard::box(id = "input_box1", title = "Peptide abundance table", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                                                  shinycssloaders::withSpinner(DT::dataTableOutput("input_body_csv_table"), type = 6, color = "#2A9C7F", size = 0.5), width = 12)),
                              shiny::fluidRow(shinydashboard::box(id = "input_box2", title = "Fasta file table", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                                                  shinycssloaders::withSpinner(DT::dataTableOutput("input_body_fasta_table"), type = 6, color = "#2A9C7F", size = 0.5), width = 12))
      ),
      #Tab 2: conditions_n_samples
      shinydashboard::tabItem(tabName = "conditions_n_samples",
                              shiny::fluidRow(shinydashboard::box(id = "cns_box1", title = "Conditions and samples table", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                                                  shinycssloaders::withSpinner(DT::dataTableOutput("cns_body_table"), type = 6, color = "#2A9C7F", size = 0.5), 
                                                                  shiny::uiOutput("cns_download_button"), width = 12))
      ),
      #Tab 3: modifications_n_normalization
      shinydashboard::tabItem(tabName = "modifications_n_normalization",
                              shiny::fluidRow(shinydashboard::box(id = "mnn_box1", title = "Modifications table", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                                                  shinycssloaders::withSpinner(DT::dataTableOutput("mnn_body_modifications_table"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                  shiny::uiOutput("mnn_modifications_download_button"),
                                                                  shiny::uiOutput("mnn_body_no_modification_text"), width = 12)),
                              shiny::fluidRow(shinydashboard::box(id = "mnn_box2", title = "Normalized data table", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                                                  shinycssloaders::withSpinner(DT::dataTableOutput("mnn_body_normalization_table"), type = 6, color = "#2A9C7F", size = 0.5) ,
                                                                  shiny::uiOutput("mnn_normalization_download_button"), width = 12))
      ),
      #Tab 4: protein_inference
      shinydashboard::tabItem(tabName = "protein_inference",
                              shiny::fluidRow(shinydashboard::box(id = "pi_box1", title = "Protein inference table", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, status = "primary",
                                                                  shinycssloaders::withSpinner(DT::dataTableOutput("pi_body_inference_table"), type = 7, color = "#2A9C7F"),
                                                                  shiny::uiOutput("pi_protein_groups_download_button"),width = 12)),
                              shiny::tags$style(shiny::HTML(".panel-title>.small, .panel-title>.small>a, .panel-title>a, .panel-title>small, .panel-title>small>a h3 {font-weight: bold;}")),
                              shiny::tags$style(shiny::HTML(".panel-title {font-size:18px;}")),
                              shiny::tags$style(shiny::HTML(".panel-primary>.panel-heading {color:#fff;background-color:#2A9C7F;} .panel-primary {border-color:#2A7F69;} .panel-primary>.panel-heading+.panel-collapse>.panel-body{border-top-color:#2A7F69;}")),
                              shinyBS::bsCollapse(id = "pi_collapse1",
                                                  shinyBS::bsCollapsePanel(title = "Connected component graph", value = "pi_panel1", style = "primary",
                                                  shiny::uiOutput("pi_body_cc_text"),
                                                  networkD3::forceNetworkOutput("pi_body_protein_graph"),
                                                  shiny::uiOutput("pi_graph_download_cui")))
      ),
      #Tab 5: protein_quantification
      shinydashboard::tabItem(tabName = "protein_quantification",
                              shiny::tags$style(".nav-tabs {background: #f4f4f4;} .nav-tabs-custom .nav-tabs li.active:hover a, .nav-tabs custom .nav-tabs li.active a {background-color: #fff; border-color: #fff;} .nav-tabs-custom .nav-tabs li.active {border-top-color: #2A7F69;}"),
                              #pq_tabbox contains all the body UI for protein_quantification tab. It is at the bottom of this R script.
                              shinycssloaders::withSpinner(shiny::uiOutput("pq_tabbox"), type = 7, color = "#2A9C7F")
      )
    )
    
  })
  ########################### END dashboardBody UI ###########################
  
  ########################### START notification dropdownMenu ###########################
  output$notification_dropdown_menu <- shinydashboard::renderMenu({

    notification_SC <- shinydashboard::notificationItem(text = "Source code", icon = icon("file-code-o"), status = "info", href = paste0("noti1"))
    notification_SC$children[[1]] <- shiny::tags$a(href = "https://bitbucket.org/vtsiamis/viqor/", target = "_blank", list(notification_SC$children[[1]]$children))
    
    notification_Man <- shinydashboard::notificationItem(text = "Manual", icon = shiny::icon("file-pdf-o"), status = "info", href = paste0("noti2"))
    notification_Man$children[[1]] <- shiny::tags$a(href = "Manual.pdf", target = "_blank", list(notification_Man$children[[1]]$children))
    
    notification_Contact <- shinydashboard::notificationItem(text = "Contact", icon = icon("envelope-square"), status = "info", href = paste0("noti3"))
    notification_Contact$children[[1]] <- shiny::tags$a(href = "mailto:vasileios@bmb.sdu.dk", target = "_blank", list(notification_Contact$children[[1]]$children))
    
    notification_Institution <- shinydashboard::notificationItem(text = "Institution", icon = icon("university"), status = "info", href = paste0("noti4"))
    notification_Institution$children[[1]] <- shiny::tags$a(href = "https://www.sdu.dk/", target = "_blank", list(notification_Institution$children[[1]]$children))
    
    notification_menu <- shinydashboard::dropdownMenu(notification_SC,
                                                      notification_Man,
                                                      notification_Contact,
                                                      notification_Institution,
                                                      type = "notifications",
                                                      badgeStatus = NULL,
                                                      icon = icon("question-circle"),
                                                      headerText = "Info:")

    return(notification_menu)
  })
  ########################### END notification dropdownMenu ###########################
  
  ########################### START INPUT TAB ###########################

  #Reactive value which indicates whether the user selected to use the example files or not.
  example_data_import <- shiny::reactiveValues(load = FALSE)
  
  #If input_example_csv_link action link is pressed then examples_data_import becomes true.
  shiny::observeEvent(input$input_example_csv_link, {
    
    example_data_import$load <- TRUE
    
  })
  
  #Reactive expression that returns the name and paths of CSV and FASTA files.
  importedData <- shiny::reactive({
    
    csv_name <- NULL
    csv_path <- NULL
    fasta_name <- NULL
    fasta_path <- NULL
    
    if(example_data_import$load){
    
      csv_name <- "Dataset/example_data.csv"
      csv_path <- "Dataset/example_data.csv"
      fasta_name <- "Dataset/example_fasta.fasta"
      fasta_path <- "Dataset/example_fasta.fasta"
      
    }
    
    if(!is.null(input$input_csv)) {
      
      csv_name <- input$input_csv$name
      csv_path <- input$input_csv$datapath
      
    }
    
    if(!is.null(input$input_fasta)) {
      
      fasta_name <- input$input_fasta$name
      fasta_path <- input$input_fasta$datapath
      
    }
    
    return(list(csv_name = csv_name, csv_path = csv_path, fasta_name = fasta_name, fasta_path = fasta_path))
    
  })
  
  #Reactive expression that loads and hold the CSV and FASTA file data. Depends on reactive expression importedData().
  inputData <- shiny::reactive({
    
    shiny::req(importedData())
    
    if(example_data_import$load){
      
      progress <- shiny::Progress$new(session, min = 0, max = 2)
      progress$set(value = 0, message = "Loading example files", detail = "CSV file...")
      
    }

    if(is.null(importedData()$csv_path)){
        
      data <- NULL
        
    } else {
        
      if(input$input_sep_csv == "tab"){
          
        error <- try(read.csv(importedData()$csv_path, header=input$input_header_csv, sep = "\t", dec = input$input_dec_csv, stringsAsFactors = FALSE, comment.char = "#"), silent = TRUE)
          
        if(class(error) != "try-error"){
          
          data <- read.csv(importedData()$csv_path, header=input$input_header_csv, sep = "\t", dec = input$input_dec_csv, stringsAsFactors = FALSE, comment.char = "#")
          
        } else {
            
          data <- NULL
            
        }
            
      } else {
          
        error <- try(read.csv(importedData()$csv_path, header = input$input_header_csv, sep = input$input_sep_csv, dec = input$input_dec_csv, stringsAsFactors = FALSE, comment.char = "#"), silent = TRUE)
          
        if(class(error) != "try-error"){
            
          data <- read.csv(importedData()$csv_path, header = input$input_header_csv, sep = input$input_sep_csv, dec = input$input_dec_csv, stringsAsFactors = FALSE, comment.char = "#")
          
        } else {
            
          data <- NULL
            
        }
          
      }
        
      if(!is.null(data)){
          
        if(ncol(data) < 2){
          
          data <- NULL
          
        }
      }
        
    }
    
    if(example_data_import$load){
      
      progress$set(value = 1, detail = "Loading fasta file...")
      
    }

    if(is.null(importedData()$fasta_path)){
        
      fasta <- NULL
      annotation <- NULL
        
    } else {
        
      fasta <- protr::readFASTA(file = importedData()$fasta_path, legacy.mode = TRUE, seqonly = FALSE)
      annotation <- seqinr::read.fasta(file = importedData()$fasta_path, seqtype = "AA", as.string = TRUE, legacy.mode = TRUE)
      annotation <- unlist(lapply(sub(" OS.*", "\\1", unlist(getAnnot(annotation))), function(x) strsplit(x, "[|]")[[1]][3]))
        
    }
      
    allData <- list(data = data, fasta = fasta, annotation = annotation)
    
    if(example_data_import$load){
      
      progress$set(value = 2, detail = "Files loaded...")
      progress$close()
      
    }
    
    return(allData)
    
  })
  
  #Conditional UI to render the name of CSV file. Initially is hidden and appears after input_go button is pressed.
  output$input_file_name_csv <- shiny::renderUI({
      
      shiny::conditionalPanel(condition = "input.sidebar_menu == 'input'",
        shiny::tags$h5(id = "input_file_names_csv_text", 
                       shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("CSV file name:")),
                       shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp;"), importedData()$csv_name)
      ))
    
  })
  
  #Conditional UI to render the name of Fasta file. Initially is hidden and appears after input_go button is pressed.
  output$input_file_name_fasta <- shiny::renderUI({
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'input'",
      shiny::tags$h5(id = "input_file_names_fasta_text",
                     shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("FASTA file name:")),
                     shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp;"), importedData()$fasta_name)
                                    
    ))
    
  })
  
  #When input$input_csv is not empty, collapse the input_box1 once.
  shiny::observeEvent(importedData()$csv_path,{
    
    js$collapse("input_box1")
    
  }, once = T)
  
  #Render CSV data table in input tab body. 
  output$input_body_csv_table <- DT::renderDataTable({
    
    shiny::req(inputData()$data)
      
    table <- inputData()$data
    
    DT::datatable(table, options = list(scrollX = T), class = "nowrap display")
    
  })
  
  #When input$input_fasta is not empty, collapse the input_box2 once.
  shiny::observeEvent(importedData()$fasta_path,{
    
    shinyjs::js$collapse("input_box2")
    
  }, once = T)
  
  #Rendering FASTA data table in input tab body. 
  output$input_body_fasta_table <- DT::renderDataTable({
    
    shiny::req(inputData()$fasta)
      
    table <- data.frame(Accession = sub(".*[|]([^.]+)[|].*", "\\1", names(inputData()$fasta)), Annotation = inputData()$annotation, Sequence = unlist(inputData()$fasta), row.names = c(1:length(unlist(inputData()$fasta))))
      
    DT::datatable(table, options = list(scrollX = T), class = "nowrap display")
    
  })
  
  #Conditional UI to render input_go action button, when there is valid input for both CSV and FASTA file.
  output$input_action <- shiny::renderUI({
    
    shiny::req(inputData()$data, inputData()$fasta)
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'input'",
                            shiny::actionButton(inputId = "input_go", label = "Next", icon = icon("fas fa-angle-right"))
    )
    
  })

  #When example_data_import$load is true, all input tab elements in sidebar hide (except input_go action button)
  shiny::observeEvent(example_data_import$load == T, {
    
    shinyjs::show(id = "input_cui_csv", anim = FALSE)
    sapply(c("input_empty1", "input_text2", "input_q_icon1", "input_csv", "input_example_csv_link", "input_q_icon1_5","input_text3", "input_text4", "input_q_icon2", "input_header_csv", "input_text5", "input_q_icon3", 
             "input_sep_csv", "input_text6", "input_q_icon4", "input_dec_csv", "input_text7", "input_q_icon5", "input_log_transformed"), function(x) shinyjs::hide(id = x, anim = FALSE))
    shinyjs::show(id = "input_cui_fasta", anim = FALSE)
    sapply(c("input_text9", "input_q_icon6", "input_fasta"), function(x) shinyjs::hide(id = x, anim = FALSE))
    
  })
  
  #When input_go action button is pressed, all the input tab elements in sidebar hide.
  #The app is then procceed to conditions_n_samples tab.  
  shiny::observeEvent(input$input_go, {
    
    #When example_data_import$load is false, all input tab elements in sidebar hide.
    if(!example_data_import$load){
    
      if(is.null(input$input_csv) & example_data_import$load == F){
        
        sapply(c("input_line1", "input_line2", "input_text1", "input_line3", "input_line4"), function(x) shinyjs::hide(id = x, anim = FALSE))
  
      } else {
        
        shinyjs::show(id = "input_cui_csv", anim = FALSE)
        
      }
      
      sapply(c("input_empty1", "input_text2", "input_q_icon1", "input_csv", "input_example_csv_link", "input_q_icon1_5","input_text3", "input_text4", "input_q_icon2", "input_header_csv", "input_text5", "input_q_icon3", 
               "input_sep_csv", "input_text6", "input_q_icon4", "input_dec_csv", "input_text7", "input_q_icon5", "input_log_transformed"), function(x) shinyjs::hide(id = x, anim = FALSE))
      
      if(is.null(input$input_fasta) & example_data_import$load == F){
        
        sapply(c("input_line5", "input_line6", "input_text8", "input_line7", "input_line8"), function(x) shinyjs::hide(id = x, anim = FALSE))
        
      } else {
        
        shinyjs::show(id = "input_cui_fasta", anim = FALSE)
        
      }
      
      sapply(c("input_text9", "input_q_icon6", "input_fasta"), function(x) shinyjs::hide(id = x, anim = FALSE))
    
    }
    
    shinyjs::hide(id = "input_go", anim = FALSE)
    
    shinydashboard::updateTabItems(session, "sidebar_menu", "conditions_n_samples")
    
  })
  ########################### END INPUT TAB ###########################
  
  ########################### START CONDITIONS AND SAMPLES TAB ###########################
  #When input_go button is pressed, the elements of conditions_n_samples tab become visible.
  shiny::observeEvent(input$input_go, {
    
    sapply(c("cns_header", "cns_blank1", "cns_tip1", "cns_quantitative_columns_slider", "cns_tip2", "cns_conditions_slider", "cns_tip3",
             "cns_reference_checkbox", "cns_tip5", "cns_arrangement_checkbox", "cns_tip4", "cns_na_slider"), function(x) shinyjs::show(id = x, anim = FALSE))
    
    if(example_data_import$load == T){
      
      value <- 6
      
    } else {
      
      value <- 4
      
    }
    
    shiny::updateSliderInput(session, inputId = "cns_total_quantitative_columns", label = NULL, min = 4, max = ncol(inputData()$data), step = 1, value = value)
    
    
  })
  
  #Adjust cns_total_conditions and cns_maximum_na sliders to the cns_total_quantitative_columns slider value.
  shiny::observeEvent(input$cns_total_quantitative_columns,{
    
    if(example_data_import$load == T){
      
      conditions.value <- 2
      maximum.na.value <- 1
      
    } else {
      
      conditions.value <- 1
      maximum.na.value <- 1
      
    }
    
    shiny::updateSliderInput(session, inputId = "cns_total_conditions", label = NULL, min = 1, max = input$cns_total_quantitative_columns, step = 1, value = conditions.value)
    shiny::updateSliderInput(session, inputId = "cns_maximum_na", label = NULL, min = 1, max = input$cns_total_quantitative_columns-2, step = 1, value = maximum.na.value)
    
  })
  
  #Conditional UI for the cns_reference_run slider when cns_q_reference checkbox is checked.
  output$cns_reference_ui <- shiny::renderUI({
    
    if(input$cns_q_reference){
      
      shiny::conditionalPanel(condition = "input.sidebar_menu == 'conditions_n_samples'",
                              shiny::sliderInput(inputId = "cns_reference_run", label = NULL, min = 1, max = 100, step = 1, value = 1, ticks = F)
      )
    }
    
  })
  
  #Adjust cns_reference_run slider maximum value to the cns_total_conditions value when cns_q_reference is true.
  shiny::observeEvent(input$cns_q_reference == T | input$cns_total_conditions,{
    
    shiny::updateSliderInput(session, inputId = "cns_reference_run", label = NULL, min = 1, max = input$cns_total_conditions, step = 1, value = 0)
    
  })
  
  #Conditional UI to render the parameters chosen by the user in conditions_n_samples tab sidebar. Initially is hidden and appears after cns_go button is pressed.
  output$cns_parameters_cui <- shiny::renderUI({
    
    shiny::req(input$cns_go)
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'conditions_n_samples'",
                            shiny::tags$h5(id = "cns_parameters_text", 
                                           shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Quantitative columns:  "), input$cns_total_quantitative_columns),
                                           shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Conditions:  "), input$cns_total_conditions),
                                           if(input$cns_q_reference == T){shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Reference condition:  "), input$cns_reference_run)},
                                           shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Maximum NA values:  "), input$cns_maximum_na)
                     ))
  })
  
  #Reactive expression that filters data based on slider values and by general.Filtering() function.
  filterData <- shiny::reactive({
    
    shiny::req(inputData()$data)
      
    if(!is.null(importedData()$fasta_path)){
       
      filtered.data <- general.Filtering(peptide.Data.Frame = inputData()$data, type = "missing", total.quantitative.columns = input$cns_total_quantitative_columns, 
                                         total.conditions = input$cns_total_conditions, maxNA = input$cns_maximum_na, log = input$input_log_transformed, column.arrangement = input$cns_arrangement)
        
    } else { #Is not used anymore since both CSV and Fasta files are required.
        
      filtered.data <- general.Filtering(peptide.Data.Frame = inputData()$data, type = "include", total.quantitative.columns = input$cns_total_quantitative_columns, 
                                         total.conditions = input$cns_total_conditions, maxNA = input$cns_maximum_na, log = input$input_log_transformed, column.arrangement = input$cns_arrangement)
      
    }
    
    return(filtered.data)
    
  })
  
  #Reactive expression that creates a structure that contains the samples grouped by condition through samples.And.Groups() function.
  labelsData <- shiny::reactive({
    
    shiny::req(filterData())
    
    if(input$cns_q_reference){
      
      shiny::req(input$cns_reference_run)    
      labels <- samples.And.Groups(peptide.Data.Frame = filterData(), total.conditions = input$cns_total_conditions, total.quantitative.columns = input$cns_total_quantitative_columns, reference.condition = input$cns_reference_run)
          
    } else {
        
      labels <- samples.And.Groups(peptide.Data.Frame = filterData(), total.conditions = input$cns_total_conditions, total.quantitative.columns = input$cns_total_quantitative_columns)
        
    }
    
    return(labels)
    
  })
  
  #When filterData() is not empty, collapse the cns_box1 once.
  shiny::observeEvent(!is.null(filterData()),{
    
    shinyjs::js$collapse("cns_box1")
    
  }, once = T)
  
  #Renders the cns_body_table in conditions_n_samples tab body.
  output$cns_body_table <- DT::renderDataTable({
    
    req(input$input_go, labelsData(), filterData())
    
    labels <- labelsData()
    filter <- filterData()
      
    if(input$cns_total_conditions >= 2){
        
      DT::datatable(filter, options = list(scrollX = T), class = "nowrap display") %>% DT::formatStyle(labels[,setdiff(seq(from =  1, to = ncol(labels), by = 2), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#FFF5FE" 
                                                                                 ) %>% DT::formatStyle(labels[,setdiff(seq(from =  2, to = ncol(labels), by = 2), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#F4FFF4"
                                                                                 ) %>% DT::formatStyle(if(input$cns_q_reference == T){labels[,input$cns_reference_run]}, backgroundColor = "#FFFFE1")
      
    } else {
        
      DT::datatable(filter, options = list(scrollX = T), class = "nowrap display") %>% DT::formatStyle(labels[,setdiff(1:ncol(labels), if(input$cns_q_reference){input$cns_reference_run == T} else {0})], backgroundColor = "#FFF5FE" 
                                                                                 ) %>% DT::formatStyle(if(input$cns_q_reference == T){labels[,input$cns_reference_run]}, backgroundColor = "#FFFFE1")
        
    }
    
  })
  
  #Download handler for the filterData() data table.
  output$cns_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("FilteredData_", Sys.time(), ".csv", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      write.csv(filterData(), file, row.names = FALSE)
      
    }
    
  )
  
  #Conditional UI to render download button in conditions_n_samples tab body when filterData() exists.
  output$cns_download_button <- shiny::renderUI({
    
    req(input$input_go, filterData())
    
    shiny::tags$p(shiny::downloadButton("cns_download", "Download"),
                  shinyBS::bsTooltip(id = "cns_download", title = "Download the above table in CSV format." , placement = "right", trigger = "hover", options = list(container = "body"))
    )

  })

  #Conditional UI to render the action button cns_go.
  output$cns_action <- shiny::renderUI({
    
    req(input$input_go)
    
    if(input$cns_total_quantitative_columns %% input$cns_total_conditions == 0){
      
      shiny::conditionalPanel(condition = "input.sidebar_menu == 'conditions_n_samples'",
                              shiny::actionButton(inputId = "cns_go", label = "Next", icon = icon("fas fa-angle-right"))
      )
      
    }
    
  })
  
  #When cns_go is pressed elements of condition_n_samples sidebar hide and application moves to modifications_n_normalization tab.
  shiny::observeEvent(input$cns_go, {
    
    shinyjs::show(id = "cns_parameters", anim = FALSE)
    sapply(c("cns_text2", "cns_q_icon1", "cns_total_quantitative_columns", "cns_text3", "cns_q_icon2", "cns_total_conditions", "cns_text4",
             "cns_q_icon3", "cns_q_reference"), function(x) shinyjs::hide(id = x, anim = FALSE))
    
    if(input$cns_q_reference){
      
      shinyjs::hide(id = "cns_reference_run", anim = FALSE)
      
    }
    
    sapply(c("cns_text6", "cns_q_icon5", "cns_arrangement", "cns_text5", "cns_q_icon4", "cns_maximum_na", "cns_go"), function(x) shinyjs::hide(id = x, anim = FALSE))
    
    shinydashboard::updateTabItems(session, "sidebar_menu", "modifications_n_normalization")
    
  })
  ########################### END CONDITIONS AND SAMPLES TAB ###########################
  
  ########################### START MODIFICATIONS AND NORMALIZATION TAB ###########################
  #When cns_go button is pressed, the elements of modifications_n_normalization tab become visible.
  shiny::observeEvent(input$cns_go, {
    
    sapply(c("mnn_header1", "mnn_blank1", "mnn_modification_types_checkbox", "mnn_header2", "mnn_tip2", "mnn_normalization_radiobuttons"), function(x) shinyjs::show(id = x, anim = FALSE))

  })
  
  #Reactive expression for the detection of modification types found in the data set by PTM.recognition().
  modificationTypes <- shiny::reactive({
    
    shiny::req(input$cns_go, filterData(), labelsData())

    types <- PTM.recognition(filterData())
      
    return(types)
    
  })
  
  #Conditional UI to render the tooltip regarding the mnn_modification_types checkbox. There are two options: modifications or no modifications.
  output$mnn_modification_types_tip <- shiny::renderUI({
    
    req(input$cns_go, modificationTypes())
    
    if(modificationTypes()$model != 0){
      
      shiny::conditionalPanel(condition = "input.sidebar_menu == 'modifications_n_normalization'",
                              shiny::tags$h5(id = "mnn_text2", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Choose modifications to separate:")),
                              shiny::tags$h6(id = "mnn_q_icon1", style = "display:inline;", icon("question-circle")),
                              shinyBS::bsTooltip(id = "mnn_q_icon1", title = "Select the modification types that should proceed to PTMs visualization analysis and press the Choose! button. The peptides of unselected modification types will be treated as unmodified peptides and therefore will be included in protein summarization process.", placement = "right", trigger = "hover", options = list(container = "body"))
      )
      
    } else {
      
      shiny::conditionalPanel(condition = "input.sidebar_menu == 'modifications_n_normalization'",
                              shiny::tags$h5(id = "mnn_text2", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("No modifications were detected.")),
                              shiny::tags$h6(id = "mnn_q_icon1", style = "display:inline;", icon("question-circle")),
                              shinyBS::bsTooltip(id = "mnn_q_icon1", title = "There are no modified peptides found in the dataset.", placement = "right", trigger = "hover", options = list(container = "body"))
      )
      
    }
    
  })
  
  #After cns_go action button is pressed and modificationTypes() is not empty, mnn_modification_types checkbox is updated to modificationTypes() values. 
  shiny::observeEvent(input$cns_go,{
    
    shiny::req(modificationTypes())
    
    if((modificationTypes()$model != 0) & (example_data_import$load == F)){
      
      shiny::updateCheckboxGroupInput(session, inputId = "mnn_modification_types", label = NULL, choices = modificationTypes()$modifications, selected = NULL)
    
    } else if((modificationTypes()$model != 0) & (example_data_import$load == T)){
      
      shiny::updateCheckboxGroupInput(session, inputId = "mnn_modification_types", label = NULL, choices = modificationTypes()$modifications, selected = modificationTypes()$modifications)
      
    } else {
      
      shinyjs::hide(id = "mnn_modification_types", anim = FALSE)
      shinyjs::show(id = "mnn_blank2", anim = FALSE)
      
    }
    
  })
  
  #Conditional UI to render action button mnn_modifications_go under mnn_modification_types checkbox when modifications exist.
  output$mnn_modifications_action <- shiny::renderUI({
    
    shiny::req(input$cns_go, modificationTypes())
    
    if(modificationTypes()$model != 0){
      
      shiny::conditionalPanel(condition = "input.sidebar_menu == 'modifications_n_normalization'",
                              shiny::actionButton(inputId = "mnn_modifications_go", label = "Choose!"))
      
    }  
    
  })
  
  #Conditional UI to render the modification types chosen by the user in modifications_n_normalization tab sidebar.
  output$mnn_parameters_cui1 <- shiny::renderUI({
    
    shiny::req(input$mnn_modifications_go)
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'modifications_n_normalization'",
                            shiny::tags$h5(id = "mnn_parameters_text1", 
                                           if(length(input$mnn_modification_types) > 1){
                                            
                                            shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp; "), shiny::tags$b("Chosen modifications:  "), 
                                            shiny::tags$p(shiny::HTML(paste(paste("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; '", input$mnn_modification_types[1:(length(input$mnn_modification_types) - 1)], "'<br>", sep = ""), collapse  = " "), paste("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; and '", input$mnn_modification_types[length(input$mnn_modification_types)], "'", sep = ""))))
                                           
                                           } else if(length(input$mnn_modification_types) == 1){
                          
                                            shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp; "), shiny::tags$b("Chosen modification:  "),
                                            shiny::tags$p(shiny::HTML(paste("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; '", input$mnn_modification_types, "'", sep = ""))))
                          
                                           } else {
                          
                                            shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp; "), shiny::tags$b("Chosen modifications:  "),
                                            shiny::tags$p(shiny::HTML(paste("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ", "none" , sep = ""))))
                          
                                           }
                        
                                         )
                          )
    
  })
  
  #Hiding of mnn_modifications_go action button and mnn_modification_types checkbox when mnn_modifications_go button is pressed.
  shiny::observeEvent(input$mnn_modifications_go,{
    
    sapply(c("mnn_text2", "mnn_q_icon1", "mnn_modification_types", "mnn_modifications_go"), function(x) shinyjs::hide(id = x, anim = FALSE))

  })
  
  #Reactive expression that uses PTM.extraction(), which creates a structure that contains all information regarding to the modifications types that have been chosen.
  modificationHandling <- shiny::reactive({
    
    req(modificationTypes(), input$mnn_modifications_go)

    modifications <- PTM.extraction(peptide.Data.Frame = filterData(), keep.modifications = input$mnn_modification_types, all.modifications = modificationTypes()$modifications, model = modificationTypes()$model)
    
    if(shiny::req(length(input$mnn_modification_types) > 0)){
      
      #Updates the colors for each modification type detected in the dataset.
      modificationsColors$colors <- data.frame(Modification = unique(unlist(modifications$type)), Color = modifications.palette[1:length(unique(unlist(modifications$type)))], stringsAsFactors = F)

    }
    
    ##### Save modifications locally #####
    #save(modifications, file="Mod.RData")
    
    return(modifications)
    
  })
  
  #Reactive value that contains the colors for each detected type of modification. Takes values through modificationHandling() reactive expression when at least 1 modification type is chosen.
  modificationsColors <- shiny::reactiveValues(colors = NA)
  
  #Reactive expression that creates the mod_table which is later rendered in mod_n_norm body.
  modificationsTable <- shiny::reactive({
    
    shiny::req(modificationHandling())
    shiny::req(length(input$mnn_modification_types) > 0)
    
    table <- data.frame(Sequence = modificationHandling()$sequences, Counterpart = modificationHandling()$counterpart)
    
    if(max(lengths(modificationHandling()$type)) > 1){
      
      modifications <- t(sapply(modificationHandling()$type, "[", seq(max(lengths(modificationHandling()$type)))))
      colnames(modifications) <- paste("Modification", 1:ncol(modifications), sep = " ")
      
      positions <- t(sapply(modificationHandling()$position, "[", seq(max(lengths(modificationHandling()$position)))))
      colnames(positions) <- paste("Position", 1:ncol(positions), sep = " ")
      
    } else {
      
      modifications <- t(t(sapply(modificationHandling()$type, "[", 1)))
      colnames(modifications) <- "Modification"
      
      positions <- t(t(sapply(modificationHandling()$position, "[", 1)))
      colnames(positions) <- "Position"
      
    }
    
    abundances <- filterData()[modificationHandling()$index, (ncol(filterData()) + 1 - input$cns_total_quantitative_columns):ncol(filterData())]
    
    table <- cbind(table, modifications, positions, abundances)
    rownames(table) <- 1:nrow(table)
    
    return(table)
    
  })
  
  #When modificationsTable() is not empty, collapse the mnn_box1 once.
  shiny::observeEvent(!is.null(modificationsTable()),{
    
    shinyjs::js$collapse("mnn_box1")
    
  }, once = T)
  
  #Text rendering in mnn_box1 instead of mnn_body_modifications_table in case where no modification was detected or in case no modification was selected.
  output$mnn_body_no_modification_text <- shiny::renderUI({
    
    
    if(shiny::req(modificationTypes()$model) == 0){
      
      shiny::tags$h4(id = "mnn_body_text1", "No modifications detected.")
      
    } else if(shiny::req(length(input$mnn_modification_types)) < 1){
      
      shiny::tags$h4(id = "mnn_body_text1", "No modification type have been selected.")
      
    }
    
  })
  
  #Rendering the modificationsTable() table in modifications_n_normalization tab body.
  output$mnn_body_modifications_table <- DT::renderDataTable({
    
    shiny::req(modificationsTable())
    
    if(input$cns_total_conditions >= 2){
      
      DT::datatable(modificationsTable(), options = list(scrollX = T), class = "nowrap display") %>% DT::formatStyle(labelsData()[,setdiff(seq(from =  1, to = ncol(labelsData()), by = 2), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#FFF5FE"
                                                                                               ) %>% DT::formatStyle(labelsData()[,setdiff(seq(from =  2, to = ncol(labelsData()), by = 2), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#F4FFF4"
                                                                                               ) %>% DT::formatStyle(if(input$cns_q_reference == T){labelsData()[,input$cns_reference_run]}, backgroundColor = "#FFFFE1")
    } else {
      
      DT::datatable(modificationsTable(), options = list(scrollX = T), class = "nowrap display") %>% DT::formatStyle(labelsData()[,setdiff(1:ncol(labelsData()), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#FFF5FE"
                                                                                               ) %>% DT::formatStyle(if(input$cns_q_reference == T){labelsData()[,input$cns_reference_run]}, backgroundColor = "#FFFFE1")
    }
    
  })
  
  #Download handler for the modifications data table.
  output$mnn_modifications_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste('ModificationsData_', Sys.time(), ".csv", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      write.csv(modificationsTable(), file, row.names = FALSE)
      
    }
    
  )
  
  #Conditional UI to render download button in conditions_n_samples tab body when modificationsTable() exists.
  output$mnn_modifications_download_button <- shiny::renderUI({
    
    shiny::req(input$mnn_modifications_go, modificationsTable())
    
    shiny::tags$p(shiny::downloadButton("mnn_modifications_download", "Download"),
                  shinyBS::bsTooltip(id = "mnn_modifications_download", title = "Download the above table in CSV format.", placement = "right", trigger = "hover", options = list(container = "body"))
    )
    
  })
  
  #Reactive value keep the filter but unnormalized unmodified data. (absolute intensities)
  unmodifiedDataAbs<- shiny::reactiveValues(data = NULL)
  
  #Reactive expression for the normalization of the unmodified peptide quantitative data.
  #If there are modifications found/selected in the dataset: remove them and normalize the remaining.
  #If there are no modifications found in the dataset: normalize the whole dataset.
  normalizedData <- shiny::reactive({
    
    shiny::req(input$cns_go, filterData(), labelsData(), modificationTypes())

    if(modificationTypes()$model != 0){
      
      shiny::req(input$mnn_modifications_go)
      
      if(length(input$mnn_modification_types) > 0){
        
        unmodified.data <- filterData()[-modificationHandling()$index,]
        rownames(unmodified.data) <- 1:nrow(unmodified.data)
        
      } else if(length(input$mnn_modification_types) == 0) {
        
        unmodified.data <- filterData()
        
      }
      
      if(modificationTypes()$model == 1){

        unmodified.data[,1] <- unlist(lapply(unmodified.data[,1], function(x) gsub("\\(.*?\\)", "", x)))
        
      } else if (modificationTypes()$model == 2){
        
        unmodified.data[,1] <- unlist(lapply(unmodified.data[,1], function(x) gsub("\\[.*?\\]", "", x)))
        
      } else if (modificationTypes()$model == 3){
        
        unmodified.data[,1] <- unlist(lapply(unmodified.data[,1], function(x) gsub("[a-z]*[a-z]", "", x)))
        
      }
      
    } else if(modificationTypes()$model == 0){
      
      unmodified.data <- filterData()
      
    } else {
      
      unmodified.data <- NULL
      
    }
    
    if(!is.null(unmodified.data)){
      
      #Residues replacement to assist protein inference. Only for unmodified peptides.
      # I <- L and J <- L replacement (J stands for either leucine or isoleucine)
      unmodified.data[,1] <- gsub("[I]|[J]", "L", unmodified.data[,1])
      
      if(input$input_log_transformed){
        
        unmodified.data[,(ncol(unmodified.data) - input$cns_total_quantitative_columns + 1):ncol(unmodified.data)] <- 2^unmodified.data[,(ncol(unmodified.data) - input$cns_total_quantitative_columns + 1):ncol(unmodified.data)]
        
      }
      
      if(length(unique(unmodified.data[,1])) != nrow(unmodified.data)){
      
        if(!is.null(input$input_fasta) | (is.null(input$input_fasta) & example_data_import$load == T)){
          
          unmodified.data <- aggregate(x = unmodified.data[,(ncol(unmodified.data) - input$cns_total_quantitative_columns + 1):ncol(unmodified.data)]/1, by = list(unmodified.data[,1]), FUN = sum, na.rm = T)
          colnames(unmodified.data)[1] <- colnames(filterData())[1]
          
        } else {
          
          unmodified.data <- aggregate(x = unmodified.data[,(ncol(unmodified.data) - input$cns_total_quantitative_columns + 1):ncol(unmodified.data)]/1, by = list(unmodified.data[,1], unmodified.data[,2]), FUN = sum, na.rm = T)
          colnames(unmodified.data)[1:2] <- colnames(filterData())[1:2]
          
        }
        
        unmodified.data[unmodified.data == 0] <- NA
      
      } else {
        
        unmodified.data <- unmodified.data[order(unmodified.data[,1]), ]
        rownames(unmodified.data) <- 1:nrow(unmodified.data)

      }
      
      #Keep the filtered but un-normalized intensities of unmodified peptides for total sum and weighted sum summarizations.
      unmodifiedDataAbs$data <- unmodified.data
      
      #Data normalization
      data <- log.Zero.Center.Normalization(peptide.Data.Frame = unmodified.data, labels = labelsData(), log = F, method = input$mnn_normalization_method)
      
    } else {
      
      data <- NULL
      
    }
    
    return(data)
    
  })
  
  #Conditional UI to render the normalization method selected by the user in modifications_n_normalization tab sidebar.
  output$mnn_parameters_cui2 <- shiny::renderUI({
    
    shiny::req(input$mnn_go)
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'modifications_n_normalization'",
                            shiny::tags$p(shiny::tags$h5(id = "mnn_text5", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Normalized by: ")),
                                          shiny::tags$h5(id = "mnn_text6", style = "display:inline;", input$mnn_normalization_method))
                        
                           )
  })
  
  #When normalizedData() is not empty, collapse the mnn_box1 once.
  shiny::observeEvent(!is.null(normalizedData()),{
    
    shinyjs::js$collapse("mnn_box2")
    
  }, once = T)
  
  #Render the normalizedData() table in modifications_n_normalization tab body.
  output$mnn_body_normalization_table <- DT::renderDataTable({
    
    shiny::req(normalizedData())
    
    if(input$cns_total_conditions >= 2){
      
      DT::datatable(normalizedData(), options = list(scrollX = T), class = "nowrap display") %>% DT::formatStyle(labelsData()[,setdiff(seq(from =  1, to = ncol(labelsData()), by = 2), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#FFF5FE"
                                                                                           ) %>% DT::formatStyle(labelsData()[,setdiff(seq(from =  2, to = ncol(labelsData()), by = 2), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#F4FFF4"
                                                                                           ) %>% DT::formatStyle(if(input$cns_q_reference == T){labelsData()[,input$cns_reference_run]}, backgroundColor = "#FFFFE1")
    } else {
      
      DT::datatable(normalizedData(), options = list(scrollX = T), class = "nowrap display") %>% DT::formatStyle(labelsData()[,setdiff(1:ncol(labelsData()), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#FFF5FE"
                                                                                           ) %>% DT::formatStyle(if(input$cns_q_reference == T){labelsData()[,input$cns_reference_run]}, backgroundColor = "#FFFFE1")
    }
    
  })
  
  #Download handler for the normalized data table.
  output$mnn_normalization_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste('NormalizedData_', input$mnn_normalization_method, "_", Sys.time(), ".csv", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      write.csv(normalizedData() , file, row.names = FALSE)
      
    }
    
  )
  
  #Conditional UI to render download button in conditions_n_samples tab body when filterData() exists.
  output$mnn_normalization_download_button <- shiny::renderUI({
    
    shiny::req(normalizedData())
    
    shiny::tags$p(shiny::downloadButton("mnn_normalization_download", "Download"),
                  shinyBS::bsTooltip(id = "mnn_normalization_download", title = "Download the above table in CSV format.", placement = "right", trigger = "hover", options = list(container = "body"))
    )
    
  })
  
  #Conditional UI for mnn_go action button render in modification_n_normalization tab sidebar, when normalizedData() exists.
  output$mnn_action <- shiny::renderUI({
    
    shiny::req(normalizedData())
      
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'modifications_n_normalization'",
                            shiny::actionButton(inputId = "mnn_go", label = "Next", icon = icon("fas fa-angle-right"))
    )
      
  })
  
  #When mnn_go action button is pressed modifications_n_normalization tab elements hide and app proceed to protein_inference tab.
  shiny::observeEvent(input$mnn_go, { 
    
    if(modificationTypes()$model == 0){
      
      shinyjs::hide("mnn_q_icon1", anim = FALSE)
      
    }
    
    sapply(c("mnn_text4", "mnn_q_icon2", "mnn_normalization_method", "mnn_go"), function(x) shinyjs::hide(id = x, anim = FALSE))    
    
    shinydashboard::updateTabItems(session, "sidebar_menu", "protein_inference")
    
  })
  ########################### END MODIFICATIONS AND NORMALIZATION TAB ###########################
  
  ########################### PROTEIN INFERENCE TAB ###########################
  #Show initial elements of protein_inference tab when mnn_go button is pressed.
  shiny::observeEvent(input$mnn_go, { 
    
    sapply(c("pi_header1", "pi_blank1", "pi_tip1", "pi_parsimony_radiobuttons"), function(x) shinyjs::show(id = x, anim = FALSE))
    
  })
  
  #When example data are selected, then the default selection for parsimony type is the strict algorithm.
  shiny::observeEvent(example_data_import$path == T, {
    
    shiny::updateRadioButtons(session, inputId = "pi_parsimony_type", label = NULL, c("soft", "strict", "no"), selected = "strict", inline = TRUE)
    
  })
  
  #Conditional UI to render pi_parsimony_go action button when pi_parsimony_type have been selected.
  output$pi_parsimony_action <- shiny::renderUI({
    
    shiny::req(input$mnn_go, input$pi_parsimony_type)
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'protein_inference'",
                            shiny::actionButton(inputId = "pi_parsimony_go", label = "Go!", icon = icon("fas fa-angle-right"))
    )
    
  })
  
  #Conditional user interface that renders the parameters selected for protein inference in protein_inference tab sidebar.
  output$pi_parameters_cui <- shiny::renderUI({
    
    shiny::req(input$pi_parsimony_go)
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'protein_inference'",
                            shiny::tags$p(shiny::HTML("&nbsp;&nbsp;&nbsp; "), shiny::tags$b("Inference by:  "), input$pi_parsimony_type)
    )
    
  })
  
  #Function to trigger elements hiding when parsimony_go action button is pressed.
  pi_parsimony_hide <- function(){
    
    sapply(c("pi_parsimony_go", "pi_parsimony_type", "pi_q_icon1", "pi_text2"), function(x) shinyjs::hide(id = x, anim = FALSE))
    
  }
  
  #Reactive expression that performs protein inference throught protein.Grouping() function when pi_parsimony_go action button is pressed.
  proteinInference <- shiny::reactive({
    
    groups <- list(proteins = NULL, peptides = NULL, cc = NULL, g = NULL, mapping = NULL)
    
    shiny::req(input$pi_parsimony_go)
    
    pi_parsimony_hide()
    
    groups <- protein.Grouping(peptide.Data.Frame =  normalizedData(), fasta = importedData()$fasta_path, parsimony = input$pi_parsimony_type)
    
    ##### Save protein groups locally #####
    #save(groups, file="Groups.RData")
    
    #If there are modified peptides that have been selected to procceed for quantification, their counterparts are mapped against the loaded FASTA file.
    if(length(input$mnn_modification_types) > 0){
      
      modifiedPeptidesMapping$map <- mapping.Modified.Peptides(counterparts = modificationHandling()$counterpart.sub, fasta = importedData()$fasta_path)

    }
    
    return(groups)
    
  })
  
  #Interactive value that keeps the modified peptides mapping.
  modifiedPeptidesMapping <- shiny::reactiveValues(map = NA)
  
  #When proteinInference()$proteins is not empty, collapse the pi_box1 once.
  shiny::observeEvent(!is.null(proteinInference()$proteins),{
    
    shinyjs::js$collapse("pi_box1")
    
  }, once = T)
  
  #Renders the pi_body_inference_table in protein_inference tab body.
  output$pi_body_inference_table <- DT::renderDataTable({
    
    if(!is.null(proteinInference()$proteins)){
      
      table <- data.frame(Protein = proteinInference()$proteins, Peptides = unlist(lapply(1:length(proteinInference()$proteins), function(x) paste(unlist(proteinInference()$peptides[x]), collapse = "|"))))
      DT::datatable(table, options = list(scrollX = T), selection = "single", class = "nowrap display")
      
    } 
    
  })
  
  #Conditional UI to render download button in protein_inference tab body when proteinInference() exists.
  output$pi_protein_groups_download_button <- shiny::renderUI({
    
    shiny::req(proteinInference())
    
    shiny::tags$p(shiny::downloadButton("pi_protein_groups_download", "Download"),
                  shinyBS::bsTooltip(id = "pi_protein_groups_download", title = "Download the above table in CSV format.", placement = "right", trigger = "hover", options = list(container = "body"))
    )
    
  })
  
  #Download handler for the pi_body_inference_table.
  output$pi_protein_groups_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("ProteinGroups_", input$pi_parsimony_type, "_", Sys.time(), ".csv", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      write.csv(data.frame(Protein = proteinInference()$proteins, Peptides = unlist(lapply(1:length(proteinInference()$proteins), function(x) paste(unlist(proteinInference()$peptides[x]), collapse = "|")))), file, row.names = FALSE)
      
    }
    
  )

  #Reactive expression that observes if any row of pi_body_inference_table is selected.
  pi_observeClick <- shiny::reactive({
    
    click <- 0
    
    if(!is.null(input$pi_body_inference_table_rows_selected)){
      
      click <- 1
      
    } else {
      
      click<- 0
      
    }
    
    return(click)
    
  })
  
  #Based on pi_observeClick observation the pi_panel1 panel collapse or expands.
  shiny:: observeEvent(pi_observeClick(), {
    
    if(pi_observeClick() == 0){
      
      shinyBS::updateCollapse(session, "pi_collapse1", close = "pi_panel1")
      
    } else {
      
      shinyBS::updateCollapse(session, "pi_collapse1", open = "pi_panel1")
      
    }
    
  })
  
  #Conditional UI for text rendering in protein_inference tab body for the connected component graphs.
  output$pi_body_cc_text <- shiny::renderUI({
    
    shiny::req(input$pi_parsimony_go)
    shiny::req(input$pi_body_inference_table_rows_selected)
    
    shiny::tags$h5(id = "pi_body_text2", paste("Connected component graph of protein group", proteinInference()$proteins[input$pi_body_inference_table_rows_selected], collapse = " "))
    
  })
  
  #Reactive expression that creates the graph of selected protein group (the whole connected component is illustrated in the graph).
  proteinGraph <- shiny::reactive({
    
    shiny::req(input$pi_parsimony_go)
    shiny::req(input$pi_body_inference_table_rows_selected)
    
    if(grepl("|", proteinInference()$proteins[input$pi_body_inference_table_rows_selected])){
      
      proteins <- unlist(strsplit(proteinInference()$proteins[input$pi_body_inference_table_rows_selected], "|", fixed = TRUE))
      
    } else {
      
      proteins <- proteinInference()$proteins[input$pi_body_inference_table_rows_selected]
      
    }
    
    peptides <- unlist(proteinInference()$peptides[input$pi_body_inference_table_rows_selected])
    
    subgraph <- igraph::induced_subgraph(proteinInference()$g, which(proteinInference()$cc$membership == proteinInference()$cc$membership[which(igraph::V(proteinInference()$g)$pname == proteins[1])]))
    
    member <- rep(0,length(igraph::V(subgraph)$pname))
    all_proteins <- igraph::V(subgraph)$pname[which(!is.na(igraph::V(subgraph)$sequence))]
    proteins_in_group <- unlist(lapply(proteins, function(x) which(igraph::V(subgraph)$pname == x)))
    member[proteins_in_group] <- "Protein in group"
    proteins_not_in_group <- unlist(lapply(setdiff(all_proteins, proteins), function(x) which(igraph::V(subgraph)$pname == x)))
    
    if(length(proteins_not_in_group) > 0){
      
      member[proteins_not_in_group] <- "Protein not in group"
      
    }
    
    all_peptides <- igraph::V(subgraph)$pname[which(is.na(igraph::V(subgraph)$sequence))]
    peptides_in_group <- unlist(lapply(peptides, function(x) which(igraph::V(subgraph)$pname == x)))
    member[peptides_in_group] <- "Peptide in group"
    peptides_not_in_group <- unlist(lapply(setdiff(all_peptides, peptides), function(x) which(igraph::V(subgraph)$pname == x)))
    
    if(length(peptides_not_in_group) > 0){
      
      member[peptides_not_in_group] <- "Peptide not in group"
      
    }
    
    subgraph.D3 <- networkD3::igraph_to_networkD3(subgraph, member)
    
    labels <- rep(NA, length(igraph::V(subgraph)$pname))
    
    labels[which(!is.na(igraph::V(subgraph)$sequence))] <- igraph::V(subgraph)$pname[which(!is.na(igraph::V(subgraph)$sequence))]
    labels[which(is.na(igraph::V(subgraph)$sequence))] <- igraph::V(subgraph)$pname[which(is.na(igraph::V(subgraph)$sequence))]
    
    size <- rep(1,length(igraph::V(subgraph)$pname))
    size[which(!is.na(igraph::V(subgraph)$sequence))] <- 10
    
    subgraph.D3$nodes[,"label"] <- labels
    subgraph.D3$nodes[,"size"] <- size
    
    plot.D3 <- networkD3::forceNetwork(Links = subgraph.D3$links, 
                                       Nodes = subgraph.D3$nodes, 
                                       Source = "source", 
                                       Target = "target", 
                                       NodeID = "label",
                                       Group = "group",
                                       Nodesize = "size",
                                       radiusCalculation = "Math.sqrt(d.nodesize)+5",
                                       colourScale = JS('d3.scaleOrdinal() .domain(["Protein in group", "Protein not in group","Peptide in group","Peptide not in group"]) .range(["#00fe00", "#fd0000","#0000ff","#ff04ff"]);'),
                                       height = 500,
                                       width = 1000,
                                       fontSize = 14,
                                       zoom = F,
                                       opacity = 1,
                                       legend = T,
                                       linkDistance = 150,
                                       bounded = T,
                                       opacityNoHover = 1,
                                       charge= -20)
    
    return(plot.D3)
    
  })
  
  #Renders the proteinGraph() in protein_inference tab body.
  output$pi_body_protein_graph <- networkD3::renderForceNetwork({
    
    shiny::req(proteinGraph())
    proteinGraph()
    
  })
  
  #Conditional UI to render pi_save_graph_go action button in protein_inference tab body. 
  output$pi_graph_download_cui <- shiny::renderUI({
    
    shiny::req(proteinInference())
    shiny::req(proteinGraph())
    shiny::req(input$pi_body_inference_table_rows_selected)

    shiny::tags$p(shiny::tags$div(style = "display: inline-block;", shiny::downloadButton("pi_graph_download_html", "Save as HTML"),
                                  shinyBS::bsTooltip(id = "pi_graph_download_html", title = "Download the above graph in HTML format.", placement = "right", trigger = "hover", options = list(container = "body"))),
                                  shiny::tags$div(style = "display: inline-block;", shiny::downloadButton("pi_graph_download_pdf", "Save as PDF"),
                                  shinyBS::bsTooltip(id = "pi_graph_download_pdf", title = "Download the above graph in PDF format.", placement = "right", trigger = "hover", options = list(container = "body"))))
    
  })
  
  #Save selected proteinGraph() in HTML format.
  output$pi_graph_download_html <- shiny::downloadHandler(
    
    filename = function() {
      
      paste(gsub("\\|", "_", proteinInference()$proteins[input$pi_body_inference_table_rows_selected]), ".html", sep = "", collapse = "")
      
    },
    content = function(file) {
      
      networkD3::saveNetwork(proteinGraph(), file, selfcontained = T)

    }
    
  )
  
  #Save selected proteinGraph() in PDF format.
  output$pi_graph_download_pdf <- shiny::downloadHandler(
    
    filename = function() {
      
      paste(gsub("\\|", "_", proteinInference()$proteins[input$pi_body_inference_table_rows_selected]), ".pdf", sep = "", collapse = "")
      
    },
    content = function(file) {
      
      temp_name <- tempfile(pattern = "ccGraph", fileext = ".html")
      
      htmlwidgets::saveWidget(proteinGraph(), temp_name)
      
      webshot2::webshot(url = temp_name, file = file)
    
      unlink(temp_name)
      unlink(paste(gsub(".html", "", temp_name), "_files", collapse = "", sep = ""), recursive = T)
      
    }
    
  )
  
  #Conditional UI to render pi_go action button to proccess in protein_quantification tab.
  output$pi_action <- shiny::renderUI({
    
    shiny::req(proteinInference())
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'protein_inference'",
                            shiny::actionButton(inputId = "pi_go", label = "Next", icon = icon("fas fa-angle-right"))
    )
    
  })
  
  #When pi_go action button is pressed protein_inference tab elements hide and appplication proceed to protein_quantification tab.
  shiny::observeEvent(input$pi_go, {
    
    shinyjs::hide(id = "pi_go", anim = FALSE)

    shinydashboard::updateTabItems(session, "sidebar_menu", "protein_quantification")
    
  })
  
  ########################### PROTEIN INFERENCE TAB END ###########################
  
  ########################### PROTEIN QUANTIFICATION TAB START ###########################
  
  #When pi_go action button is pressed protein_quantification sidebar elements appear.
  shiny::observeEvent(input$pi_go, {
    
    sapply(c("pq_header1", "pq_blank1", "pq_tip1", "pq_farms_weight_slider", "pq_farms_mu_slider"), function(x) shinyjs::show(id = x, anim = FALSE))
    
  })
  
  #When cns_q_reference is checked and a reference run is chosen then this conditional UI is rerdered in the protein_quantification tab to give both options of rescale probes (reference, average).
  output$pq_probes_rescale_cui <- shiny::renderUI({
    
    shiny::req(input$pi_go, input$cns_q_reference)
    
    if(input$cns_reference_run != 0){
      shiny::conditionalPanel(condition = "input.sidebar_menu == 'protein_quantification'",
                              shiny::tags$p(shiny::tags$h5(id = "pq_text3", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Relative quantification: ")),
                                            shiny::tags$h6(id = "pq_q_icon2", style = "display:inline;", icon("question-circle"))),
                                            shinyBS::bsTooltip(id = "pq_q_icon2", title = "Select the relative quantification method that should be applied on the peptide abundances of each protein group.", placement = "right", trigger = "hover", options = list(container = "body")),
                                            shiny::radioButtons(inputId = "pq_farms_rescale", label = NULL, choices = c("By average", "By reference"), selected = "By average", inline = TRUE)
      )
    }
    
  })
  
  #Conditional UI to generate pq_quantification_go action button.
  output$pq_quantification_action <- shiny::renderUI({
    
    shiny::req(input$pi_go)
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'protein_quantification'",
                            shiny::actionButton(inputId = "pq_quantification_go", label = "Go!", icon = icon("fas fa-angle-right"))
    )
    
  })
  
  #Conditional UI to render the fast-FARMS parameters chosen by the user in protein_inference sidebar.
  output$pq_parameters_cui <- shiny::renderUI({
    
    shiny::req(input$pq_quantification_go)
    
    shiny::conditionalPanel(condition = "input.sidebar_menu == 'protein_quantification'",
                            shiny::tags$p(shiny::tags$h5(id = "pq_text4", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("fast-FARMS weight: ")),
                                          shiny::tags$h5(id = "pq_text5", style = "display:inline;", input$pq_farms_weight)),
                            shiny::tags$p(shiny::tags$h5(id = "pq_text6", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("fast-FARMS mu: ")),
                                          shiny::tags$h5(id = "pq_text7", style = "display:inline;", input$pq_farms_mu)),
                            shiny::tags$p(shiny::tags$h5(id = "pq_text8", style = "display:inline;", shiny::HTML("&nbsp;&nbsp;&nbsp;"), shiny::tags$b("Rescale: "),
                                          if(input$cns_q_reference == T){shiny::tags$h5(id = "pq_text9", style = "display:inline;", input$pq_farms_rescale)} else {shiny::tags$h5(id = "pq_text9", style = "display:inline;", "By average")}))
    )
    
  })
  
  #When pq_quantification_go action button is pressed sliders and input elements hide from the protein_quantification sidebar by the following function.
  pq_quantification_hide <- function(){
    
    
    sapply(c("pq_text2", "pq_q_icon1"), function(x) shinyjs::hide(id = x, anim = FALSE))
    
    if(input$cns_q_reference){
      if(input$cns_reference_run != 0){
        
        sapply(c("pq_text3", "pq_q_icon2", "pq_farms_rescale"), function(x) shinyjs::hide(id = x, anim = FALSE))
        
      }
    }
    
    sapply(c("pq_farms_weight", "pq_farms_mu", "pq_text4", "pq_q_icon3", "pq_quantification_go"), function(x) shinyjs::hide(id = x, anim = FALSE))
    
  }
  
  #Reactive expression to calculate all the weights for peptide probes cooresponding to all protein groups.
  probesAndWeights <- shiny::reactive({
    
    probes_n_weights <- list(protein.group = NULL, probes = NULL, weights = NULL)
    
    shiny::req(input$pq_quantification_go)
    
    pq_quantification_hide()
    
    cat("\nProtein quantification started.\n")
    cat("| Weights calculation started.\n")
    cat("| |\n")
    
    if(input$cns_q_reference){
      
      if(input$cns_reference_run != 0){
        
        shiny::req(input$pq_farms_rescale)
        
        if(input$pq_farms_rescale == "By average"){
          
          probes_n_weights <- probes.And.Weights(peptide.abundances = normalizedData(), protein.groups =  proteinInference(), labels =  labelsData(), rescale = T, weight = input$pq_farms_weight, mu = input$pq_farms_mu)
          
        } else if(input$pq_farms_rescale == "By reference"){
          
          probes_n_weights <- probes.And.Weights(peptide.abundances =  normalizedData(), protein.groups =  proteinInference(), labels =  labelsData(), rescale = F, weight = input$pq_farms_weight, mu = input$pq_farms_mu)
          
        }
        
      }
      
    } else {
      
      probes_n_weights <- probes.And.Weights(peptide.abundances =  normalizedData(), protein.groups = proteinInference(), labels =  labelsData(), rescale = T, weight = input$pq_farms_weight, mu = input$pq_farms_mu)
      
    }
    
    cat("| Weights calculated.\n")
    
    ##### Save probes and weights locally #####
    #save(probes_n_weights, file="PnW.RData")
    
    return(probes_n_weights)
    
  })
  
  #Reactive expression to calculate the protein expression for all protein groups for different fast-FARMS weight threshold cutoffs. 
  quantificationThreshold <- shiny::reactive({
    
    shiny::req(input$pq_quantification_go, probesAndWeights())
    
    cat("|\n")
    cat("| Protein expression calculation started.\n")
    
    shiny::withProgress(message = "Protein quantification: ", min = 0, max = (length(proteinInference()$proteins) + 1), value = length(proteinInference()$proteins), detail = "It is going to take some time...", {
      
      shiny::incProgress(amount = length(proteinInference()$proteins))
      
      expr <- calculateParallel(probes_n_weights =  probesAndWeights(), cores = 11)
      
    })

    cat("| |\n")
    cat("| Protein expression calculated.\n")
    cat("Protein quantification finished.\n")
    
    ##### Save quantificationThreshold() locally #####
    #save(expr, file="Expr.RData")
    
    return(expr)
    
  })
  
  #Slider in the body of protein_quantification tab. It is used to navigate through different weight cutoffs.
  output$pq_body_weight_slider <- shiny::renderUI({
    
    shiny::req(input$pq_quantification_go, quantificationThreshold())
    
    shiny::tags$div(style = "margin: 0px 20px 0px 0px;", shiny::sliderInput(inputId = "pq_body_weight", label = "Choose the peptide weight threshold:", min = 0, max = 1, step = 0.1, value = 0, ticks = F))
    
  })
  
  #Reactive expression that calculates the protein quantitative data of the selected weight threshold per condition.
  quantificationCondition <- shiny::reactive({
    
    shiny::req(input$pq_quantification_go, quantificationThreshold())

    index <- input$pq_body_weight*10 + 1
    expr <- sapply(c(1:ncol(labelsData())), function(x) rowMeans(quantificationThreshold()[[index]][,labelsData()[,x]], na.rm = T))
      
    expr <- data.frame(probesAndWeights()$protein.group, expr, stringsAsFactors = F)
    colnames(expr) <- c("Protein", colnames(labelsData()))
    
    return(expr)
    
  })
  
  #Reactive expression to calculate the protein expression by weighted sum approach. 
  quantificationThresholdAbs <- shiny::reactive({
    
    shiny::req(input$pq_quantification_go, probesAndWeights(), unmodifiedDataAbs$data)

    expr <- quantification.By.Summation(peptide.abundances = unmodifiedDataAbs$data, labels = labelsData(), log = F, normalization.method = input$mnn_normalization_method, probes_n_weights = probesAndWeights())

    return(expr)
    
  })
  
  #Reactive value to map the quantified proteins to their annotations retrieved from reference fasta file.
  proteinAnnotation <- shiny::reactiveValues(annotation = NULL)
  
  #Protein annotation is filled once.
  shiny::observeEvent(c(input$pq_quantification_go, quantificationThreshold()), {
    
    first.protein.ids <- sub("\\|.*", "\\1", quantificationThreshold()[[1]][,1])
    progress <- shiny::Progress$new(session, min = 0, max = length(first.protein.ids))
    progress$set(value = 0, message = "Protein annotations retrieve:")
    annotation <- inputData()$annotation[unlist(sapply(1:length(first.protein.ids), function(x){ which = which(sub(".*[|]([^.]+)[|].*", "\\1", names(inputData()$fasta)) == first.protein.ids[x])
                                                                                                 progress$set(value = x, detail = paste("Annotation: ", x, " of ", length(first.protein.ids), sep = "", collapse = ""))
                                                                                                 return(which)}))]
    progress$close()
    proteinAnnotation$annotation <- annotation
    
  })
  
  #Reactive expression to count how many peptides are considered for protein summarization, for every protein and for every possible value of input$pq_body_weight
  peptideCounter <- shiny::reactive({
    
    shiny::req(probesAndWeights(), quantificationThreshold())
    
    if(input$cns_total_quantitative_columns != input$cns_total_conditions){
      
      shiny::req(quantificationCondition())
      
    }
    
    counts <- lapply(seq(0,1,0.1), function(y) unlist(lapply(1:length(probesAndWeights()$protein.group), function(x) sum(unlist(probesAndWeights()$weights[x]) >= y))))
    
    return(counts)
    
  })
  
  #Calculate and save correlation tables for all weight threshold expression tables and then creates the global correlation indicator(s).
  #Returns the graphs and the order of smiley face icons.
  globalCorrelationIndicator <- shiny::reactive({
    
    shiny::req(quantificationThreshold())
    
    expr <- quantificationThreshold()
    
    cor.data <- lapply(1:11, function(x) cor(expr[[x]][,2:ncol(expr[[x]])], use = "pairwise.complete.obs", method = "pearson"))
    
    #Save correlation matrices for all weight cut off values locally
    #save(cor.data, file="cor.RData")

    if(nrow(labelsData()) == 1){
      
      same.condition <- rep(NA, 11)
      dif.condition <- unlist(lapply(cor.data, function(x) mean(abs(x[lower.tri(x)]))))
      smiley.order <- NULL
      rect.max.min <- c(max(dif.condition, na.rm = TRUE), min(dif.condition, na.rm = TRUE))
      
    } else {
      
      boolean.matrix <- matrix(FALSE, nrow = input$cns_total_quantitative_columns, ncol = input$cns_total_quantitative_columns)

      for (i in seq(from = 1, to = input$cns_total_quantitative_columns, by = nrow(labelsData()))) {
        
        boolean.matrix[i:(i+(nrow(labelsData())-1)), i:(i+(nrow(labelsData())-1))] <- TRUE
        
      }
      
      mapping.matrix <- lower.tri(cor.data[[1]], diag = FALSE) - boolean.matrix
      diag(mapping.matrix) <- -2
      which.same <- which(mapping.matrix == -1)
      which.different <- which(mapping.matrix == 1)

      same.condition <- unlist(lapply(1:11, function(x) mean(cor.data[[x]][which.same], na.rm = TRUE)))
      dif.condition <- unlist(lapply(1:11, function(x) mean(abs(cor.data[[x]][which.different]), na.rm = TRUE)))
      order <- order(same.condition, decreasing = TRUE)
      smiley.order <- rep(NA,11)
      
      for (i in 1:11) {
        
        smiley.order[order[i]] <- c(1,2,2,3,3,4,4,5,5,6,6)[i]
        
      }
      
      rect.max.min <- c(max(same.condition, na.rm = TRUE), min(same.condition, na.rm = TRUE))
      
    }
    
    m <- list(
      l = 10,
      r = 30,
      b = 0,
      t = 0,
      pad = 0
    )
    
    p <- plotly::plot_ly(width = 280, height = 80) %>% 
                         plotly::layout(autosize = F, margin = m, showlegend = FALSE, xaxis = list(fixedrange = TRUE, tick0 = 0, dtick = 0.1, tickmode = "linear", tickfont = list(size = 9)), 
                                                                                      yaxis = list(fixedrange = TRUE, tickfont = list(size = 9)),
                                                                                      yaxis2 = list(fixedrange = TRUE, tickfont = list(size = 9), side = "right", overlaying = "y"))
    
    if(all(!is.na(same.condition))){
      
      p <- plotly::add_trace(p, x = seq(from = 0, to = 1, by = 0.1), y = same.condition, name = "Within",type = "scatter", mode = "lines", yaxis = "y")
      p <- plotly::add_trace(p, x = seq(from = 0, to = 1,by = 0.1), y = dif.condition, name = "Between", type = "scatter", mode = "lines", yaxis = "y2", opacity = 0.5)
      
    } else {
      
      p <- plotly::add_trace(p, x = seq(from = 0, to = 1,by = 0.1), y = dif.condition, name = "Between", type = "scatter", mode = "lines", yaxis = "y")
      
    }
    
    p <- plotly::config(p, displayModeBar = FALSE)
    
    plist <- vector("list", 11)
    
    for (i in 1:11) {
      
      if(i == 1){
        
        x.zero <- (i-1)/10 
        x.one <- (i-1)/10 + 0.01
        
      } else if(i == 11){
        
        x.zero <- (i-1)/10 - 0.01
        x.one <- (i-1)/10
        
      } else {
        
        x.zero <- (i-1)/10 - 0.01
        x.one <- (i-1)/10 + 0.01
        
      }
      
      plist[[i]] <- p %>% layout(shapes = list(list(type = "rect",
                                                    fillcolor = "green", line = list(color = "green"), opacity = 0.3,
                                                    x0 = x.zero, x1 = x.one, xref = "x",
                                                    y0 = rect.max.min[2], 
                                                    y1 = rect.max.min[1], yref = "y")))
      
    }
    
    return(list(graph = plist, order = smiley.order))
                      
  })
  
  #Renders the correlation indicator to protein quantification body.
  output$pq_correlation_indicator <- plotly::renderPlotly({
    
    shiny::req(globalCorrelationIndicator(), input$pq_body_weight)
    
    return(globalCorrelationIndicator()$graph[[input$pq_body_weight*10 + 1]])
    
  })
  
  #Renders the corresponding smiley icon next to correlation indicator plot according to the global correlation value.
  output$pq_smileyImage <- shiny::renderImage({
    
    shiny::req(globalCorrelationIndicator()$order, input$pq_body_weight)
    
    file <- paste("smileys/", globalCorrelationIndicator()$order[input$pq_body_weight*10 + 1], ".svg", collapse = "", sep = "")
    
    return(list(src = file, contentType = "image/svg+xml", alt = globalCorrelationIndicator()$order[input$pq_body_weight*10 + 1]))

  }, deleteFile = FALSE)
  
  ########################### VISUALIZATION MODULES FOR THRESHOLD TAB ###########################
  
  ###################################### TABLE ######################################
  ###################################################################################
  
  #Reactive value that keeps the state of the threshold expr datatable object when pq_body_weight slider is adjusted to different weight threshold.
  pq_threshold_table_state <- shiny::reactiveValues(order = list(), selected = NA, page = NA, start = 0, length = 10, search = "")
  
  #Render quantificationThreshold() protein expression table for selected weight cutoff value. 
  output$pq_body_threshold_table <- DT::renderDataTable({
    
    shiny::req(probesAndWeights(), quantificationThreshold(), input$pq_body_weight, peptideCounter(), proteinAnnotation$annotation)

    table <- quantificationThreshold()[[(input$pq_body_weight * 10 + 1)]]
    #Addition of total peptides for each group and the peptides that were considered for summarization.
    table$'Total identified peptides' <- peptideCounter()[[1]]
    table$'Total quantified peptides' <- peptideCounter()[[(input$pq_body_weight * 10 + 1)]]
    #Addition of annotation for each protein.
    table$"Annotation" <- proteinAnnotation$annotation

    table <- table[, c(1, 2:4 + input$cns_total_quantitative_columns, c(2:(input$cns_total_quantitative_columns + 1)))]
    
    #Addition of signal to noise ratio per identified protein group.
    table$"S/N" <- round(10 * log10((1 - probesAndWeights()$noise)/probesAndWeights()$noise), digits = 2)
    
    #Create color gradient for S/N column.
    breaks <- quantile(table$"S/N", probs = seq(.05, .95, .05), na.rm = TRUE)
    color.gradient <- colorRampPalette(c("#F4ADAD", "#FDF7CD", "#CEF0A1"))
    
    if(input$cns_total_conditions >= 2){
      
      DT::datatable(table, options = list(scrollX = T, order = pq_threshold_table_state$order, stateSave = TRUE, stateDuration = -1, displayStart = pq_threshold_table_state$start, pageLength = pq_threshold_table_state$length, search = list(`search` = pq_threshold_table_state$search)),
                    selection = list(mode = "single", selected = pq_threshold_table_state$selected), class = "nowrap display") %>% DT::formatStyle(labelsData()[,setdiff(seq(from =  1, to = ncol(labelsData()), by = 2), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#FFF5FE"
                                                                                                                             ) %>% DT::formatStyle(labelsData()[,setdiff(seq(from =  2, to = ncol(labelsData()), by = 2), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#F4FFF4"
                                                                                                                             ) %>% DT::formatStyle(if(input$cns_q_reference == T){labelsData()[,input$cns_reference_run]}, backgroundColor = "#FFFFE1"
                                                                                                                             ) %>% DT::formatStyle("S/N", backgroundColor = styleInterval(breaks, color.gradient(20)))
    } else {
      
      DT::datatable(table, options = list(scrollX = T, order = pq_threshold_table_state$order, stateSave = TRUE, stateDuration = -1, displayStart = pq_threshold_table_state$start, pageLength = pq_threshold_table_state$length),
                    selection = list(mode = "single", selected = pq_threshold_table_state$selected), class = "nowrap display") %>% DT::formatStyle(labelsData()[,setdiff(1:ncol(labelsData()), if(input$cns_q_reference == T){input$cns_reference_run})], backgroundColor = "#FFF5FE"
                                                                                                                             ) %>% DT::formatStyle(if(input$cns_q_reference == T){labelsData()[,input$cns_reference_run]}, backgroundColor = "#FFFFE1"
                                                                                                                             ) %>% DT::formatStyle("S/N", backgroundColor = styleInterval(breaks, color.gradient(20)))
    }
    
  })
  
  #When pq_body_weight slider changes value of weight threshold, pq_threshold_table_state values change to load selected threshold table on that state.
  shiny::observeEvent(input$pq_body_weight, {
    
    pq_threshold_table_state$order <- input$pq_body_threshold_table_state$order
    pq_threshold_table_state$selected <- input$pq_body_threshold_table_rows_selected
    pq_threshold_table_state$page <- (input$pq_body_threshold_table_state$start/input$pq_body_threshold_table_state$length + 1)
    pq_threshold_table_state$start <- input$pq_body_threshold_table_state$start
    pq_threshold_table_state$length <- input$pq_body_threshold_table_state$length
    pq_threshold_table_state$search <- if(is.null(input$pq_body_threshold_table_state$`search`$`search`)){""} else {input$pq_body_threshold_table_state$`search`$`search`}

  })
  
  #Conditional UI to render download button in protein quantification threshold tab when quantificationThreshold() is not empty.
  output$pq_threshold_table_download_button <- shiny::renderUI({
    
    shiny::req(quantificationThreshold())
    
    shiny::tags$p(shiny::downloadButton("pq_threshold_table_download", "Download"),
                  shinyBS::bsTooltip(id = "pq_threshold_table_download", title = "Download the above table in CSV format.", placement = "right", trigger = "hover", options = list(container = "body"))
    )
    
  })
  
  #Download handler for the protein expression data table in threshold tab.
  output$pq_threshold_table_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("ProteinExpression_", input$pq_body_weight, "_",Sys.time(), ".csv", colapse = "", sep = "")
      
    },
    content = function(file) {
      
      shiny::req(input$pq_body_weight)
      
      table <- quantificationThreshold()[[(input$pq_body_weight * 10 + 1)]]
      table$'Total identified peptides' <- peptideCounter()[[1]]
      table$'Total quantified peptides' <- peptideCounter()[[(input$pq_body_weight * 10 + 1)]]
      table$"Annotation" <- proteinAnnotation$annotation
      table <- table[, c(1, 2:4 + input$cns_total_quantitative_columns, c(2:(input$cns_total_quantitative_columns + 1)))]
      table$"S/N" <- round(10 * log10((1 - probesAndWeights()$noise)/probesAndWeights()$noise), digits = 2)
      
      write.csv(table, file, row.names=FALSE)
      
    }
    
  )
  
  #Conditional UI to render download button for weighted summation approach in protein quantification threshold tab.
  output$pq_weighted_sum_sample_download_button <- shiny::renderUI({
    
    shiny::req(quantificationThreshold(), quantificationThresholdAbs())
    
    shiny::tags$p(shiny::downloadButton("pq_weighted_sum_sample_download", "Download (Weighted Sum)"),
                  shinyBS::bsTooltip(id = "pq_weighted_sum_sample_download", title = "Download the protein expressions by weighted summation summarization in CSV format.", placement = "left", trigger = "hover", options = list(container = "body"))
    )
    
  })
  
  #Download handler for the protein expression data table summarized by weighted sum in threshold tab.
  output$pq_weighted_sum_sample_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("ProteinExpressionWeightedSum_", Sys.time(), ".csv", colapse = "", sep = "")
      
    },
    content = function(file) {
      
      shiny::req(quantificationThresholdAbs())
      
      table <- quantificationThresholdAbs()$weighted
      
      write.csv(table, file, row.names=FALSE)
      
    }
    
  )
  
  #Conditional UI to render download button for total summation summarization in protein quantification threshold tab.
  output$pq_total_sum_sample_download_button <- shiny::renderUI({
    
    shiny::req(quantificationThreshold(), quantificationThresholdAbs())
    
    shiny::tags$p(shiny::downloadButton("pq_total_sum_sample_download", "Download (Total Sum)"),
                  shinyBS::bsTooltip(id = "pq_total_sum_sample_download", title = "Download the protein expressions by total summation summarization in CSV format.", placement = "left", trigger = "hover", options = list(container = "body"))
    )
    
  })
  
  #Download handler for the protein expression data table summarized by total sum in threshold tab.
  output$pq_total_sum_sample_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("ProteinExpressionTotalSum_", Sys.time(), ".csv", colapse = "", sep = "")
      
    },
    content = function(file) {
      
      shiny::req(quantificationThresholdAbs())
      
      table <- quantificationThresholdAbs()$total
      
      write.csv(table, file, row.names=FALSE)
      
    }
    
  )
  
  #################################### LINE PLOT ####################################
  ###################################################################################
  
  #Reactive expression that observes if any row of pq_body_threshold_table is selected.
  pq_threshold_observeClick <- shiny::reactive({
    
    click <- 0
    
    if(!is.null(input$pq_body_threshold_table_rows_selected)){
      
      click <- 1
      
    } else {
      
      click<- 0
      
    }
    
    return(click)
    
  })
  
  #Based on pq_threshold_observeClick observation the pq_panel1 panel collapse or expands.
  shiny::observeEvent(pq_threshold_observeClick(), {
    
    if(pq_threshold_observeClick() == 0){
      
      shinyBS::updateCollapse(session, "pq_collapse1", close = "pq_panel1")
      
    } else {
      
      shinyBS::updateCollapse(session, "pq_collapse1", open = "pq_panel1")
      
    }
    
  })
  
  #Reactive expression that creates a line plot for a chosen protein group of pq_body_threshold_table.
  lineplotThreshold <- shiny::reactive({
    
    shiny::req(quantificationThreshold(), input$pq_body_threshold_table_rows_selected, input$pq_body_weight)

    index <- input$pq_body_threshold_table_rows_selected

    p <- interactive.Line.Plot(probes =  probesAndWeights()$probes[[index]], weights = probesAndWeights()$weights[[index]], expr = quantificationThreshold()[[(input$pq_body_weight*10 + 1)]][index,], threshold = input$pq_body_weight, colors = palette, xlabel = "Runs", peptide_mapping = proteinInference()$mapping)
    
    return(p)
    
  })
  
  #Renders the interactive line plot of selected protein group in threshold tab.
  output$pq_threshold_lineplot <- plotly::renderPlotly({
    
    shiny::req(quantificationThreshold(), lineplotThreshold())
    
    return(lineplotThreshold())
    
  })
  
  #Conditional UI to create download button for the interactive line plot of the protein selected in threshold tab in a static version.
  output$pq_threshold_lineplot_download_cui <- shiny::renderUI({
    
    shiny::req(quantificationThreshold(), input$pq_body_threshold_table_rows_selected, lineplotThreshold())

    shiny::tags$p(shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_threshold_lineplot_width", label = "Width: ", min = 100, max = 1000, value = 1000, step = 100, width = "250px")),
                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_threshold_lineplot_height", label = "Height: ", min = 100, max = 1000, value = 500, step = 100, width = "200px")),
                  shiny::tags$div(style = "display: inline-block;margin: 25px 0px 0px 0px;", shiny::downloadButton("pq_threshold_lineplot_download", "Download")),
                  shinyBS::bsTooltip(id = "pq_threshold_lineplot_download", title = "Choose the width and height and download the above figure in PDF format.", placement = "right", trigger = "hover", options = list(container = "body")))
    
  })
  
  #Download handler that creates and downloads in PDF format the interactive line plot of the protein selected in threshold tab.
  output$pq_threshold_lineplot_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("LinePlot", gsub("\\|","_",probesAndWeights()$protein.group[input$pq_body_threshold_table_rows_selected]), ".pdf", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      temp_name <- tempfile(pattern = "lineplot", fileext = ".html")
      
      p <- lineplotThreshold()
      
      if(input$pq_threshold_lineplot_width > 0 & input$pq_threshold_lineplot_height > 0){
        
        p$width <- input$pq_threshold_lineplot_width
        p$height <- input$pq_threshold_lineplot_height
        
        for (i in 3:length(p$x$attrs)) {
          
          p$x$attrs[[i]]$opacity <- 1
          
        }
        
        htmlwidgets::saveWidget(p, temp_name)
        
        webshot2::webshot(url = temp_name, file = file)
        
        unlink(temp_name)
        unlink(paste(gsub(".html", "", temp_name), "_files", collapse = "", sep = ""), recursive = T)
        
        shiny::updateNumericInput(session, "pq_threshold_lineplot_width", value = 1000)
        shiny::updateNumericInput(session, "pq_threshold_lineplot_height", value = 500)
        
      } 
      
    }
    
  )
  
  ######################################### VIQoR PLOT #########################################
  ###############################################################################################
  
  #Based on pq_threshold_observeClick observation the pq_panel3 panel collapse or expands.
  shiny::observeEvent(pq_threshold_observeClick(), {
    
    if(pq_threshold_observeClick() == 0){
      
      shinyBS::updateCollapse(session, "pq_collapse3", close = "pq_panel3")
      
    } else {
      
      shinyBS::updateCollapse(session, "pq_collapse3", open = "pq_panel3")
      
    }
    
  })
  
  #Selectinput objects pq_threshold_select1 and pq_threshold_select2 choices are updated to the non-NA column names of the row selected in pq_body_threshold_table.
  shiny::observeEvent(input$pq_body_threshold_table_rows_selected, {
    
    table <- quantificationThreshold()[[(input$pq_body_weight * 10 + 1)]]
    na_columns <- which(is.na(table[input$pq_body_threshold_table_rows_selected,]))
    
    if(!all(is.na(table[input$pq_body_threshold_table_rows_selected, 2:ncol(table)]))){
    
      if(length(na_columns) > 0){
      
        table <- table[, -which(is.na(table[input$pq_body_threshold_table_rows_selected, ]))]
      
      }
    
      shiny::updateSelectInput(session, "pq_threshold_select1", choices = colnames(table)[2:ncol(table)], selected = colnames(table)[2])
      shiny::updateSelectInput(session, "pq_threshold_select2", choices = colnames(table)[2:ncol(table)], selected = colnames(table)[2])
    
    } else {
      
      shiny::updateSelectInput(session, "pq_threshold_select1", choices = NULL, selected = NULL)
      shiny::updateSelectInput(session, "pq_threshold_select2", choices = NULL, selected = NULL)
      
    }
    
  })
  
  #Reactive expression that creates the VIQoR plot for a chosen protein group and two comparing runs.
  VIQoRPlotThreshold <- shiny::reactive({
   
    shiny::req(input$pq_body_weight, probesAndWeights(), quantificationThreshold(), proteinInference(), input$pq_body_threshold_table_rows_selected, input$pq_threshold_select1, input$pq_threshold_select2, input$pq_threshold_enzyme)
      
    if(length(modifiedPeptidesMapping$map) == 1){
      
      p <- VIQoR.Plot(protein_group_index =  input$pq_body_threshold_table_rows_selected, 
                      condition1 = input$pq_threshold_select1,
                      condition2 = input$pq_threshold_select2,
                      protein_expression = quantificationThreshold(),
                      probes_n_weights = probesAndWeights(),
                      threshold = input$pq_body_weight,
                      peptide_mapping = proteinInference()$mapping,
                      pep_colors =  palette,
                      enzyme = input$pq_threshold_enzyme)
      
    } else {
        
      p <- VIQoR.Plot(protein_group_index = input$pq_body_threshold_table_rows_selected,
                      condition1 = input$pq_threshold_select1, 
                      condition2 = input$pq_threshold_select2, 
                      protein_expression = quantificationThreshold(), 
                      probes_n_weights = probesAndWeights(), 
                      threshold = input$pq_body_weight, 
                      peptide_mapping = proteinInference()$mapping, 
                      pep_colors =  palette, 
                      modified_peptide_mapping = modifiedPeptidesMapping$map, 
                      modifications = modificationHandling(), 
                      data = if(input$input_log_transformed){filterData()[,-c(1:(ncol(filterData())-input$cns_total_quantitative_columns))]} else {log2(filterData()[,-c(1:(ncol(filterData())-input$cns_total_quantitative_columns))])},
                      mod_colors =  modificationsColors$colors, 
                      enzyme = input$pq_threshold_enzyme)  
        
    }
      
    return(p)
      
  })
  
  #Visualize the VIQoR plot in the threshold tab.
  output$pq_threshold_VIQoR <- plotly::renderPlotly({
  
    shiny::req(VIQoRPlotThreshold())
    
    return(VIQoRPlotThreshold())
  
  })
  
  #Button to download the interactive VIQoR plot of the protein selected in threshold tab in a static version.
  output$pq_threshold_VIQoR_download_cui <- shiny::renderUI({
    
    shiny::req(quantificationThreshold(), input$pq_body_threshold_table_rows_selected, VIQoRPlotThreshold())

    shiny::tags$p(shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_threshold_VIQoR_width", label = "Width: ", min = 100, max = 1000, value = 1000, step = 100, width = "250px")),
                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_threshold_VIQoR_height", label = "Height: ", min = 100, max = 1000, value = 500, step = 100, width = "200px")),
                  shiny::tags$div(style = "display: inline-block;margin: 25px 0px 0px 0px;", shiny::downloadButton("pq_threshold_VIQoR_download", "Download")),
                  shinyBS::bsTooltip(id = "pq_threshold_VIQoR_download", title = "Choose the width and height and download the above figure in PDF format.", placement = "right", trigger = "hover", options = list(container = "body")))

  })
  
  #Save selected VIQoR plot of threshold panel in PDF format.
  output$pq_threshold_VIQoR_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("VIQoRPlot_", gsub("\\|", "_", probesAndWeights()$protein.group[input$pq_body_threshold_table_rows_selected]), "_", input$pq_threshold_select1, "_vs_", input$pq_threshold_select2, ".pdf", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      temp_name <- tempfile(pattern = "VIQoR", fileext = ".html")
      
      p <- VIQoRPlotThreshold()
      
      if(input$pq_threshold_VIQoR_height > 0 & input$pq_threshold_VIQoR_width > 0){
        
        p$width <- input$pq_threshold_VIQoR_width
        p$height <- input$pq_threshold_VIQoR_height
        
        for (i in 3:length(p$x$attrs)) {
          
          p$x$attrs[[i]]$opacity <- 1
          
        }
      
        htmlwidgets::saveWidget(p, temp_name)
        
        webshot2::webshot(url = temp_name, file = file, vwidth = input$pq_threshold_VIQoR_width, vheight = input$pq_threshold_VIQoR_height, cliprect = "viewport")
        
        unlink(temp_name)
        unlink(paste(gsub( ".html", "", temp_name), "_files", collapse = "", sep = ""), recursive = T)
        
        shiny::updateNumericInput(session, "pq_threshold_VIQoR_width", value = 1000)
        shiny::updateNumericInput(session, "pq_threshold_VIQoR_height", value = 500)
        
      }
      
    }
    
  )
  
  ########################### PROTEIN EXPRESSION HEATMAP ############################
  ###################################################################################
  
  #Updates the choices of pq_threshold_n_clusters selectizeInput object.
  shiny::observeEvent(input$pq_quantification_go,{
    
    if(input$cns_total_quantitative_columns <= 2){
      
      shiny::updateSelectizeInput(session, "pq_threshold_n_clusters", label = "Select number of clusters:", choices = c(1:input$cns_total_quantitative_columns), selected = 1)
      
    } else {
      
      shiny::updateSelectizeInput(session, "pq_threshold_n_clusters", label = "Select number of clusters:", choices = c("Any", 1:input$cns_total_quantitative_columns), selected = "Any")
      
    }
    
  })
  
  #Reactive expression that creates the protein expression correlation heatmap in threshold tab.
  proteinHeatmapThreshold <- shiny::reactive({
    
    shiny::req(quantificationThreshold(), input$pq_threshold_clust_select, input$pq_threshold_distance_select, input$pq_threshold_n_clusters, input$pq_threshold_gradient, input$pq_body_weight)

    p <- correlation.Heatmap(expr =  quantificationThreshold()[[input$pq_body_weight * 10 + 1]],
                             clust_method = input$pq_threshold_clust_select, 
                             distance_method = input$pq_threshold_distance_select,
                             n_cluster = input$pq_threshold_n_clusters,
                             gradient = input$pq_threshold_gradient,
                             threshold = input$pq_body_weight, 
                             label = "Run labels")
    
    return(p)
    
  })
  
  #Renders proteinHeatmapThreshold() to threshold tab.
  output$pq_threshold_protein_heatmap <- plotly::renderPlotly({
    
    shiny::req(proteinHeatmapThreshold())
    
    return(proteinHeatmapThreshold())
    
  })
  
  #Render the download button for protein expression correlation heatmap in static version.
  output$pq_threshold_protein_heatmap_download_cui <- shiny::renderUI({
    
    shiny::req(quantificationThreshold(), proteinHeatmapThreshold())
    
    shiny::tags$p(shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_threshold_protein_heatmap_width", label = "Width: ", min = 100, max = 1000, value = 1000, step = 100, width = "250px")),
                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_threshold_protein_heatmap_height", label = "Height: ", min = 100, max = 1000, value = 500, step = 100, width = "200px")),
                  shiny::tags$div(style = "display: inline-block;margin: 25px 0px 0px 0px;", shiny::downloadButton("pq_threshold_protein_heatmap_download", "Download")),
                  shinyBS::bsTooltip(id = "pq_threshold_protein_heatmap_download", title = "Choose the width and height and download the above figure in PDF format.", placement = "right", trigger = "hover", options = list(container = "body")))
    
  })
  
  #Save protein expression correlation heatmap of threshold tab in PDF format.
  output$pq_threshold_protein_heatmap_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("ProteinExpressionHeatmap_", input$pq_body_weight, "_",Sys.time(), ".pdf", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      temp_name <- tempfile(pattern = "proteinHeatmap", fileext = ".html")
      
      p <- proteinHeatmapThreshold()
      
      if(input$pq_threshold_protein_heatmap_width > 0 & input$pq_threshold_protein_heatmap_height > 0){
        
        p$width <- input$pq_threshold_protein_heatmap_width
        p$height <- input$pq_threshold_protein_heatmap_height
        
        htmlwidgets::saveWidget(p, temp_name)
        
        webshot2::webshot(url = temp_name, file = file)
        
        unlink(temp_name)
        unlink(paste(gsub(".html", "", temp_name), "_files", collapse = "", sep = ""), recursive = T)
        
        shiny::updateNumericInput(session, "pq_threshold_protein_heatmap_width", value = 1000)
        shiny::updateNumericInput(session, "pq_threshold_protein_heatmap_height", value = 500)
        
      }
      
    }
    
  )
  
  ########################### PEPTIDE ABUNDANCE HEATMAP #############################
  ###################################################################################

  #Based on pq_threshold_observeClick observation the pq_panel5 panel collapse or expands.
  shiny::observeEvent(pq_threshold_observeClick(), {
    
    if(pq_threshold_observeClick() == 0){
      
      shinyBS::updateCollapse(session, "pq_collapse5", close = "pq_panel5")
      
    } else {
      
      shinyBS::updateCollapse(session, "pq_collapse5", open = "pq_panel5")
      
    }
    
  })
  
  #Updates the choices of pq_threshold_n_clusters_pep selectizeInput object.
  shiny::observeEvent(input$pq_body_threshold_table_rows_selected,{
    
    if(peptideCounter()[[input$pq_body_weight * 10 + 1]][input$pq_body_threshold_table_rows_selected] <= 2){

      shiny::updateSelectizeInput(session, "pq_threshold_n_clusters_pep", label = "Select number of clusters:", choices = c(1:peptideCounter()[[input$pq_body_weight * 10 + 1]][input$pq_body_threshold_table_rows_selected]), selected = 1)


    } else {

      shiny::updateSelectizeInput(session, "pq_threshold_n_clusters_pep", label = "Select number of clusters:", choices = c("Any", 1:peptideCounter()[[input$pq_body_weight * 10 + 1]][input$pq_body_threshold_table_rows_selected]), selected = "Any")

    }
    
  })

  #Reactive expression that creates the peptide expression correlation heatmap in threshold tab.
  peptideHeatmapThreshold <- shiny::reactive({
    
    shiny::req(probesAndWeights(), quantificationThreshold(), input$pq_body_threshold_table_rows_selected, input$pq_body_weight, input$pq_threshold_n_clusters_pep)

    probes <- probesAndWeights()$probes[[input$pq_body_threshold_table_rows_selected]]
    probes <- probes[which(probesAndWeights()$weights[[input$pq_body_threshold_table_rows_selected]] >= input$pq_body_weight), ]
    expr <- t(probes[, 2:ncol(probes)])
    colnames(expr) <- probes[, 1] 

    p <- correlation.Heatmap(expr = expr, 
                            clust_method = input$pq_threshold_clust_select_pep,
                            distance_method = input$pq_threshold_distance_select_pep,
                            n_cluster = input$pq_threshold_n_clusters_pep,
                            gradient = input$pq_threshold_gradient_pep,
                            threshold = input$pq_body_weight,
                            label = "Peptides")
    return(p)
    
  })
  
  #Conditional UI when the peptide correlation heatmap cannot be generated.
  output$pq_threshold_no_peptide_heatmap <- shiny::renderUI({
    
    shiny::req(input$pq_body_threshold_table_rows_selected, is.null(peptideHeatmapThreshold()))
    
    shiny::tags$p(shiny::tags$h5(id = "pq_threshold_no_peptide_heatmap_text", style = "display:inline", "Peptide correlation heatmap cannot be generated for this protein group."),
                  shiny::tags$h6(id = "pq_threshold_no_peptide_heatmap_icon", style = "display:inline", icon("question-circle")),
                  shinyBS::bsTooltip(id = "pq_threshold_no_peptide_heatmap_icon", title = "Too many missing values in peptide abundance vectors, which successively lead to correlation vectors with less than 2 finite values or a single peptide remained to that protein group.", placement = "right", trigger = "hover", options = list(container = "body")))
    
  })
  
  #Renders peptideHeatmapThreshold() to threshold tab.
  output$pq_threshold_peptide_heatmap <- plotly::renderPlotly({
    
    shiny::req(peptideHeatmapThreshold())
    return(peptideHeatmapThreshold())

  })
  
  #Renders download button for peptide abundance correlation heatmap in static version.
  output$pq_threshold_peptide_heatmap_download_cui <- shiny::renderUI({
    
    shiny::req(quantificationThreshold(), probesAndWeights(), peptideHeatmapThreshold())
    
    if(!is.null(peptideHeatmapThreshold())){
    
      shiny::tags$p(shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_threshold_peptide_heatmap_width", label = "Width: ", min = 100, max = 1000, value = 1000, step = 100, width = "250px")),
                    shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_threshold_peptide_heatmap_height", label = "Height: ", min = 100, max = 1000, value = 500, step = 100, width = "200px")),
                    shiny::tags$div(style = "display: inline-block;margin: 25px 0px 0px 0px;", shiny::downloadButton("pq_threshold_peptide_heatmap_download", "Download")),
                    shinyBS::bsTooltip(id = "pq_threshold_peptide_heatmap_download", title = "Choose the width and height and download the above figure in PDF format.", placement = "right", trigger = "hover", options = list(container = "body")))
    }
      
  })
  
  #Save peptide abundance correlation heatmap in PDF format.
  output$pq_threshold_peptide_heatmap_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("PeptideAbundanceHeatmap_", gsub("\\|", "_", probesAndWeights()$protein.group[input$pq_body_threshold_table_rows_selected]), "_",Sys.time(), ".pdf", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      temp_name <- tempfile(pattern = "peptideHeatmap", fileext = ".html")
      
      p <- peptideHeatmapThreshold()
      
      if(input$pq_threshold_peptide_heatmap_width > 0 & input$pq_threshold_peptide_heatmap_height > 0){
        
        p$width <- input$pq_threshold_peptide_heatmap_width
        p$height <- input$pq_threshold_peptide_heatmap_height
        
        htmlwidgets::saveWidget(p, temp_name)
        
        webshot2::webshot(url = temp_name, file = file)
        
        unlink(temp_name)
        unlink(paste(gsub(".html", "", temp_name), "_files", collapse = "", sep = ""), recursive = T)
        
        shiny::updateNumericInput(session, "pq_threshold_peptide_heatmap_width", value = 1000)
        shiny::updateNumericInput(session, "pq_threshold_peptide_heatmap_height", value = 500)
        
      }
      
    }
    
  )
  
  ########################### VISUALIZATION MODULES FOR THRESHOLD TAB END ###########################
  
  ########################### VISUALIZATION MODULES FOR CONDITIONS TAB ###########################
  
  ###################################### TABLE ######################################
  ###################################################################################
  
  #Reactive value that keeps the state of the condition expr datatable object in condition tab when pq_body_weight slider is adjusted to different weight threshold.
  pq_condition_table_state <- shiny::reactiveValues(order = list(), selected = NA, page = NA, start = 0, length = 10, search = "")
  
  #Renders quantificationCondition() protein expression table for the weight selected in pq_body_weight slider.
  output$pq_body_condition_table <- DT::renderDataTable({
    
    shiny::req(quantificationCondition(), peptideCounter(), proteinAnnotation$annotation)
    
    table <- quantificationCondition()
    table$'Total identified peptides' <- peptideCounter()[[1]]
    table$'Total quantified peptides' <- peptideCounter()[[input$pq_body_weight * 10 + 1]]
    #Addition of annotation for each protein.
    table$"Annotation" <- proteinAnnotation$annotation
    
    table <- table[, c(1, 2:4 + input$cns_total_conditions, 2:(input$cns_total_conditions + 1))]
    
    #Addition of signal to noise ratio per identified protein group.
    table$"S/N" <- round(10 * log10((1 - probesAndWeights()$noise)/probesAndWeights()$noise), digits = 2)
    
    #Create color gradient for S/N column.
    breaks <- quantile(table$"S/N", probs = seq(.05, .95, .05), na.rm = TRUE)
    color.gradient <- colorRampPalette(c("#F4ADAD", "#FDF7CD", "#CEF0A1"))
    
    DT::datatable(table, options = list(scrollX = T, order = pq_condition_table_state$order, stateSave = TRUE, stateDuration = -1, displayStart = pq_condition_table_state$start, pageLength = pq_condition_table_state$length, search = list(search = pq_condition_table_state$search)),
                  selection = list(mode = "single", selected = pq_condition_table_state$selected), class = "nowrap display") %>% DT::formatStyle(setdiff(colnames(labelsData())[seq(from =  1, to = ncol(labelsData()), by = 2)], if(input$cns_q_reference == T){colnames(labelsData())[input$cns_reference_run]}), backgroundColor = "#FFF5FE"
                                                                                                                           ) %>% DT::formatStyle(setdiff(colnames(labelsData())[seq(from =  2, to = ncol(labelsData()), by = 2)], if(input$cns_q_reference == T){colnames(labelsData())[input$cns_reference_run]}), backgroundColor = "#F4FFF4"
                                                                                                                           ) %>% DT::formatStyle(if(input$cns_q_reference == T){colnames(labelsData())[input$cns_reference_run]}, backgroundColor = "#FFFFE1"
                                                                                                                           ) %>% DT::formatStyle("S/N", backgroundColor = styleInterval(breaks, color.gradient(20)))
    
  })
  
  #When pq_body_weight slider changes value of weight threshold, pq_condition_table_state values change to load selected threshold table on that state.
  shiny::observeEvent(input$pq_body_weight, {
    
    pq_condition_table_state$order <- input$pq_body_condition_table_state$order
    pq_condition_table_state$selected <- input$pq_body_condition_table_rows_selected
    pq_condition_table_state$page <- (input$pq_body_condition_table_state$start/input$pq_body_condition_table_state$length + 1)
    pq_condition_table_state$start <- input$pq_body_condition_table_state$start
    pq_condition_table_state$length <- input$pq_body_condition_table_state$length
    pq_condition_table_state$search <- if(is.null(input$pq_body_condition_table_state$`search`$`search`)){""} else {input$pq_body_condition_table_state$`search`$`search`}
    
  })
  
  #Conditional UI to render the download button in protein quantification condition tab when quantificationThreshold() is not empty.
  output$pq_condition_table_download_button <- shiny::renderUI({
    
    shiny::req(quantificationCondition())
    
    shiny::tags$p(shiny::downloadButton("pq_condition_table_download", "Download"),
                  shinyBS::bsTooltip(id = "pq_condition_table_download", title = "Download the above table in CSV format.", placement = "right", trigger = "hover", options = list(container = "body"))
    )
    
  })
  
  #Download handler for the protein expression data table in condition tab.
  output$pq_condition_table_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste('ProteinExpressionConditions_', input$pq_body_weight, "_",Sys.time(), '.csv', collapse = "", sep = "")
      
    },
    content = function(file) {
      
      table <- quantificationCondition()
      table$'Total identified peptides' <- peptideCounter()[[1]]
      table$'Total quantified peptides' <- peptideCounter()[[input$pq_body_weight * 10 + 1]]
      table$"Annotation" <- proteinAnnotation$annotation
      table <- table[, c(1, 2:4 + input$cns_total_conditions, 2:(input$cns_total_conditions + 1))]
      table$"S/N" <- round(10 * log10((1 - probesAndWeights()$noise)/probesAndWeights()$noise), digits = 2)
      
      write.csv(table, file, row.names=FALSE)
      
    }
    
  )
  
  #Conditional UI to render the download button for the weighted summation summarization per condition.
  output$pq_weighted_sum_condition_download_button <- shiny::renderUI({
    
    shiny::req(quantificationThreshold(), quantificationThresholdAbs())
    
    shiny::tags$p(shiny::downloadButton("pq_weighted_sum_condition_download", "Download (Weighted Sum)"),
                  shinyBS::bsTooltip(id = "pq_weighted_sum_condition_download", title = "Download the protein expressions by weighted summation summarization in CSV format.", placement = "left", trigger = "hover", options = list(container = "body"))
    )
    
  })
  
  #Download handler for the pq_weighted_sum_condition_download_button button.
  output$pq_weighted_sum_condition_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste('ProteinExpressionWeightedSumConditions_', Sys.time(), '.csv', collapse = "", sep = "")
      
    },
    content = function(file) {
      
      shiny::req(quantificationThresholdAbs())
      
      table <- sapply(c(1:ncol(labelsData())), function(x) rowMeans(quantificationThresholdAbs()$weighted[,labelsData()[,x]], na.rm = T))
      table <- data.frame(quantificationThresholdAbs()$weighted[,1], table, stringsAsFactors = F)
      colnames(table) <- c("Protein", colnames(labelsData()))
      
      write.csv(table, file, row.names=FALSE)
      
    }
    
  )
  
  #Conditional UI to render the download button for the total summation summarization per condition.
  output$pq_total_sum_condition_download_button <- shiny::renderUI({
    
    shiny::req(quantificationThreshold(), quantificationThresholdAbs())
    
    shiny::tags$p(shiny::downloadButton("pq_total_sum_condition_download", "Download (Total Sum)"),
                  shinyBS::bsTooltip(id = "pq_total_sum_condition_download", title = "Download the protein expressions by total summation summarization in CSV format.", placement = "left", trigger = "hover", options = list(container = "body"))
    )
    
  })
  
  #Download handler for the pq_total_sum_condition_download_button button.
  output$pq_total_sum_condition_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste('ProteinExpressionTotalSumConditions_', Sys.time(), '.csv', collapse = "", sep = "")
      
    },
    content = function(file) {
      
      shiny::req(quantificationThresholdAbs())
      
      table <- sapply(c(1:ncol(labelsData())), function(x) rowMeans(quantificationThresholdAbs()$total[,labelsData()[,x]], na.rm = T))
      table <- data.frame(quantificationThresholdAbs()$total[,1], table, stringsAsFactors = F)
      colnames(table) <- c("Protein", colnames(labelsData()))
      
      write.csv(table, file, row.names=FALSE)
      
    }
    
  )

  #################################### LINE PLOT ####################################
  ###################################################################################
  
  #Reactive expression that observes if any row of pq_body_condition_table is selected.
  pq_condition_observeClick <- shiny::reactive({
    
    click <- 0
    
    if(!is.null(input$pq_body_condition_table_rows_selected)){
      
      click <- 1
      
    } else {
      
      click<- 0
      
    }
    
    return(click)
    
  })
  
  #Based on pq_threshold_observeClick observation the pq_panel1 panel collapse or expands.
  shiny::observeEvent(pq_condition_observeClick(), {
    
    if(pq_condition_observeClick() == 0){
      
      shinyBS::updateCollapse(session, "pq_collapse2", close = "pq_panel2")
      
    } else {
      
      shinyBS::updateCollapse(session, "pq_collapse2", open = "pq_panel2")
      
    }
    
  })
  
  #Reactive expression that creates the line plot for a chosen protein group of pq_body_condition_table in condition tab.
  lineplotCondition <- shiny::reactive({
    
    shiny::req(quantificationCondition(), input$pq_body_weight, input$pq_body_condition_table_rows_selected)
    
    threshold <- input$pq_body_weight
    
    index <- input$pq_body_condition_table_rows_selected
    
    probes <- probesAndWeights()$probes[[index]]
    condition.probes <- sapply(c(1:ncol(labelsData())), function(x) rowMeans(probes[, labelsData()[, x]], na.rm = T))
    colnames(condition.probes) <- colnames(quantificationCondition()[index, 2:ncol(quantificationCondition())])
    probes <- data.frame(probes[, 1], condition.probes)
    
    p <- interactive.Line.Plot(probes = probes, weights = probesAndWeights()$weights[[index]], expr = quantificationCondition()[index,], threshold = threshold, colors = palette,  xlabel = "Conditions", peptide_mapping = proteinInference()$mapping)
    
    return(p)
    
  })
  
  #Renders the interactive line plot of selected protein group in condition tab.
  output$pq_condition_lineplot <- plotly::renderPlotly({
    
    shiny::req(quantificationCondition(), lineplotCondition())
    
    return(lineplotCondition())
    
  })
  
  #Button to download the interactive line plot of the protein selected in pq_body_condition_table in a static version.
  output$pq_condition_lineplot_download_cui <- shiny::renderUI({
    
    shiny::req(quantificationCondition(), input$pq_body_condition_table_rows_selected, lineplotCondition())
    
    shiny::tags$p(shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_condition_lineplot_width", label = "Width: ", min = 100, max = 1000, value = 1000, step = 100, width = "250px")),
                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_condition_lineplot_height", label = "Height: ", min = 100, max = 1000, value = 500, step = 100, width = "200px")),
                  shiny::tags$div(style = "display: inline-block;margin: 25px 0px 0px 0px;", shiny::downloadButton("pq_condition_lineplot_download", "Download")),
                  shinyBS::bsTooltip(id = "pq_condition_lineplot_download", title = "Choose the width and height and download the above figure in PDF format.", placement = "right", trigger = "hover", options = list(container = "body")))
    
  })
  
  #Download handler that creates and downloads in PDF format the interactive line plot of the protein selected in pq_body_condition_table.
  output$pq_condition_lineplot_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("LinePlotCondition_", gsub("\\|","_",probesAndWeights()$protein.group[input$pq_body_condition_table_rows_selected]), ".pdf", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      temp_name <- tempfile(pattern = "lineplot", fileext = ".html")
      
      p <- lineplotCondition()
      
      if(input$pq_condition_lineplot_width > 0 & input$pq_condition_lineplot_height > 0){
        
        p$width <- input$pq_condition_lineplot_width
        p$height <- input$pq_condition_lineplot_height
        
        for (i in 3:length(p$x$attrs)) {
          
          p$x$attrs[[i]]$opacity <- 1
          
        }
        
        htmlwidgets::saveWidget(p, temp_name)
        
        webshot2::webshot(url = temp_name, file = file)
        
        unlink(temp_name)
        unlink(paste(gsub(".html", "", temp_name), "_files", collapse = "", sep = ""), recursive = T)
        
        shiny::updateNumericInput(session, "pq_condition_lineplot_width", value = 1000)
        shiny::updateNumericInput(session, "pq_condition_lineplot_height", value = 500)
        
      }
      
    }
    
  )
  
  ######################################### VIQoR PLOT ##########################################
  ###############################################################################################
  
  #Based on pq_condition_observeClick observation the pq_panel4 panel collapse or expands.
  shiny::observeEvent(pq_condition_observeClick(), {
    
    if(pq_condition_observeClick() == 0){
      
      shinyBS::updateCollapse(session, "pq_collapse4", close = "pq_panel4")
      
    } else {
      
      shinyBS::updateCollapse(session, "pq_collapse4", open = "pq_panel4")
      
    }
    
  })
  
  #Selectinput objects pq_condition_select1 and pq_condition_select2 choices are updated to the non-NA column names of the row selected in pq_body_condition_table.
  shiny::observeEvent(input$pq_body_condition_table_rows_selected, {
    
    table <- quantificationCondition()
    na_columns <- which(is.na(table[input$pq_body_condition_table_rows_selected,]))
    
    if(!all(is.na(table[input$pq_body_condition_table_rows_selected, 2:ncol(table)]))){
      
      if(length(na_columns) > 0){
        
        table <- table[, -which(is.na(table[input$pq_body_condition_table_rows_selected, ]))]
        
      }
      
      shiny::updateSelectInput(session, "pq_condition_select1", choices = colnames(table)[2:ncol(table)], selected = colnames(table)[2])
      shiny::updateSelectInput(session, "pq_condition_select2", choices = colnames(table)[2:ncol(table)], selected = colnames(table)[2])
      
    } else {
      
      shiny::updateSelectInput(session, "pq_condition_select1", choices = NULL, selected = NULL)
      shiny::updateSelectInput(session, "pq_condition_select2", choices = NULL, selected = NULL)
      
    }
    
  })
  
  #Reactive expression that creates the VIQoR plot for a chosen protein group and two comparing conditions.
  VIQoRPlotCondition <- shiny::reactive({
    
    shiny::req(input$pq_body_weight, probesAndWeights(), quantificationCondition(), proteinInference(), input$pq_body_condition_table_rows_selected, input$pq_condition_select1, input$pq_condition_select2)
    
    if(length(modifiedPeptidesMapping$map) == 1){
        
      p <- VIQoR.Plot(protein_group_index = input$pq_body_condition_table_rows_selected, 
                      condition1 = input$pq_condition_select1, 
                      condition2 = input$pq_condition_select2, 
                      protein_expression = quantificationCondition(), 
                      probes_n_weights = probesAndWeights(), 
                      threshold = input$pq_body_weight, 
                      peptide_mapping = proteinInference()$mapping, 
                      pep_colors = palette,
                      labels = labelsData(), 
                      enzyme = input$pq_condition_enzyme)
        
    } else {
        
      expr <- data.frame(sapply(c(1:ncol(labelsData())), function(x) rowMeans(filterData()[, labelsData()[, x]], na.rm = T)), stringsAsFactors = F)
      colnames(expr) <- colnames(labelsData())
        
      p <- VIQoR.Plot(protein_group_index = input$pq_body_condition_table_rows_selected, 
                      condition1 = input$pq_condition_select1, 
                      condition2 = input$pq_condition_select2, 
                      protein_expression = quantificationCondition(), 
                      probes_n_weights = probesAndWeights(), 
                      threshold = input$pq_body_weight, 
                      peptide_mapping = proteinInference()$mapping, 
                      pep_colors = palette, 
                      modified_peptide_mapping = modifiedPeptidesMapping$map, 
                      modifications = modificationHandling(), 
                      data = if(input$input_log_transformed){expr} else {log2(expr)},
                      mod_colors = modificationsColors$colors, 
                      labels = labelsData(), 
                      enzyme = input$pq_condition_enzyme)  
        
    }
      
    return(p)
      
  })
  
  #Visualize the VIQoR plot in the gantt plot panel in condition tab.
  output$pq_condition_VIQoR <- plotly::renderPlotly({
    
    shiny::req(VIQoRPlotCondition())
    
    return(VIQoRPlotCondition())
    
  })
  
  #Button to download the interactive VIQoR plot of the protein selected in condition tab in a static version
  output$pq_condition_VIQoR_download_cui <- shiny::renderUI({
    
    shiny::req(quantificationCondition(), input$pq_body_condition_table_rows_selected, VIQoRPlotCondition())
    
    shiny::tags$p(shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_condition_VIQoR_width", label = "Width: ", min = 100, max = 1000, value = 1000, step = 100, width = "250px")),
                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_condition_VIQoR_height", label = "Height: ", min = 100, max = 1000, value = 500, step = 100, width = "200px")),
                  shiny::tags$div(style = "display: inline-block;margin: 25px 0px 0px 0px;", shiny::downloadButton("pq_condition_VIQoR_download","Download")),
                  shinyBS::bsTooltip(id = "pq_condition_VIQoR_download", title = "Choose the width and height and download the above figure in PDF format.", placement = "right", trigger = "hover", options = list(container = "body")))
    
  })
  
  #Save selected VIQoR plot of condition panel in PDF format.
  output$pq_condition_VIQoR_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("VIQoRPlot_", gsub("\\|", "_",probesAndWeights()$protein.group[input$pq_body_condition_table_rows_selected]), "_", input$pq_condition_select1, "_vs_", input$pq_condition_select2,".pdf", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      temp_name <- tempfile(pattern = "VIQoR", fileext = ".html")
      
      p <- VIQoRPlotCondition()
      
      if(input$pq_condition_VIQoR_width > 0 & input$pq_condition_VIQoR_height > 0){
        
        p$width <- input$pq_condition_VIQoR_width
        p$height <- input$pq_condition_VIQoR_height
        
        for (i in 3:length(p$x$attrs)) {
          
          p$x$attrs[[i]]$opacity <- 1
          
        }
        
        htmlwidgets::saveWidget(p, temp_name)
        
        webshot2::webshot(url = temp_name, file = file)
        
        unlink(temp_name)
        unlink(paste(gsub( ".html", "", temp_name), "_files", collapse = "", sep = ""), recursive = T)
        
        shiny::updateNumericInput(session, "pq_condition_gantt_width", value = 1000)
        shiny::updateNumericInput(session, "pq_condition_gantt_height", value = 500)
        
      }
      
    }
    
  )
  
  ########################### PROTEIN EXPRESSION HEATMAP ############################
  ###################################################################################
  
  #Updates the choices of pq_condition_n_clusters selectizeInput object.
  shiny::observeEvent(input$pq_quantification_go,{

    if(input$cns_total_conditions <= 2){
      
      shiny::updateSelectizeInput(session, "pq_condition_n_clusters", label = "Select number of clusters:", choices = c(1:input$cns_total_conditions), selected = 1)
      
    } else {
      
      shiny::updateSelectizeInput(session, "pq_condition_n_clusters", label = "Select number of clusters:", choices = c("Any",1:input$cns_total_conditions), selected = "Any")
      
    }
    
  })
  
  #Reactive expression that creates the protein expression correlation heatmap in correlation tab.
  proteinHeatmapCondition <- shiny::reactive({
    
    shiny::req(quantificationCondition(), input$pq_condition_clust_select, input$pq_condition_distance_select, input$pq_condition_n_clusters, input$pq_condition_gradient, input$pq_body_weight)

    p <- correlation.Heatmap(expr =  quantificationCondition(),
                             clust_method = input$pq_condition_clust_select,
                             distance_method =  input$pq_condition_distance_select,
                             n_cluster = input$pq_condition_n_clusters,
                             gradient = input$pq_condition_gradient,
                             threshold = input$pq_body_weight,
                             label = "Condition labels")
    
    return(p)
    
  })
  
  #Renders proteinHeatmapCondition() to condition tab.
  output$pq_condition_protein_heatmap <- plotly::renderPlotly({
    
    shiny::req(proteinHeatmapCondition())
    
    return(proteinHeatmapCondition())
    
  })
  
  #Renders the download button for protein expression correlation heatmap of condition tab in static version.
  output$pq_condition_protein_heatmap_download_cui <- shiny::renderUI({
    
    shiny::req(quantificationCondition(), proteinHeatmapCondition())
    
    shiny::tags$p(shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_condition_protein_heatmap_width", label = "Width: ", min = 100, max = 1000, value = 1000, step = 100, width = "250px")),
                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_condition_protein_heatmap_height", label = "Height: ", min = 100, max = 1000, value = 500, step = 100, width = "200px")),
                  shiny::tags$div(style = "display: inline-block;margin: 25px 0px 0px 0px;", downloadButton("pq_condition_protein_heatmap_download", "Download")),
                  shinyBS::bsTooltip(id = "pq_condition_protein_heatmap_download", title = "Choose the width and height and download the above figure in PDF format.", placement = "right", trigger = "hover", options = list(container = "body")))
    
  })
  
  #Save protein expression correlation heatmap of condition tab in PDF format.
  output$pq_condition_protein_heatmap_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("ConditionProteinExpressionHeatmap_", Sys.time(), ".pdf", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      temp_name <- tempfile(pattern = "proteinHeatmap", fileext = ".html")
      
      p <- proteinHeatmapCondition()
      
      if(input$pq_condition_protein_heatmap_width > 0 & input$pq_condition_protein_heatmap_height > 0){
        
        p$width <- input$pq_condition_protein_heatmap_width
        p$height <- input$pq_condition_protein_heatmap_height
        
        htmlwidgets::saveWidget(p, temp_name)
        
        webshot2::webshot(url = temp_name, file = file)
        
        unlink(temp_name)
        unlink(paste(gsub(".html", "", temp_name), "_files", collapse = "", sep = ""), recursive = T)
        
        shiny::updateNumericInput(session, "pq_condition_protein_heatmap_width", value = 1000)
        shiny::updateNumericInput(session, "pq_condition_protein_heatmap_height", value = 500)
        
      }
      
    }
    
  )
  
  ########################### PEPTIDE ABUNDANCE HEATMAP #############################
  ###################################################################################
  
  #Based on pq_condition_observeClick observation the pq_panel6 panel collapse or expands.
  shiny::observeEvent(pq_condition_observeClick(), {
    
    if(pq_condition_observeClick() == 0){
      
      shinyBS::updateCollapse(session, "pq_collapse6", close = "pq_panel6")
      
    } else {
      
      shinyBS::updateCollapse(session, "pq_collapse6", open = "pq_panel6")
      
    }
    
  })
  
  #Updates the choices of pq_condition_n_clusters_pep selectizeInput object.
  shiny::observeEvent(input$pq_body_condition_table_rows_selected,{
    
    if(peptideCounter()[[input$pq_body_weight * 10 + 1]][input$pq_body_condition_table_rows_selected] <= 2){
      
      shiny::updateSelectizeInput(session, "pq_condition_n_clusters_pep", label = "Select number of clusters:", choices = c(1:peptideCounter()[[input$pq_body_weight * 10 + 1]][input$pq_body_condition_table_rows_selected]), selected = 1)
      
    } else {
      
      shiny::updateSelectizeInput(session, "pq_condition_n_clusters_pep", label = "Select number of clusters:", choices = c("Any", 1:peptideCounter()[[input$pq_body_weight * 10 + 1]][input$pq_body_condition_table_rows_selected]), selected = "Any")
      
    }
    
  })
  
  #Reactive expression that creates the peptide abundance correlation heatmap in condition tab.
  peptideHeatmapCondition <- shiny::reactive({
    
    shiny::req(probesAndWeights(), quantificationCondition(), input$pq_body_condition_table_rows_selected, input$pq_body_weight)

    probes <- probesAndWeights()$probes[[input$pq_body_condition_table_rows_selected]]
    probes <- probes[which(probesAndWeights()$weights[[input$pq_body_condition_table_rows_selected]] >= input$pq_body_weight), ]
    
    peptides <- probes[, 1]
    probes <- sapply(c(1:ncol(labelsData())), function(x) rowMeans(probes[, labelsData()[, x]], na.rm = T))
    
    if(length(peptides == 1)){
      
      expr <- data.frame(pep <- t(probes))
      
    } else {
      
      expr <- t(probes)
      
    }

    if(ncol(expr) > 1){

      colnames(expr) <- peptides
      p <- correlation.Heatmap(expr = expr, clust_method = input$pq_condition_clust_select_pep, distance_method = input$pq_condition_distance_select_pep, n_cluster = input$pq_condition_n_clusters_pep, gradient = input$pq_condition_gradient_pep, threshold = input$pq_body_weight, label = "Peptides")
    
    } else {
      
      p <- NULL
      
    }
      
    return(p)
    
  })
  
  #Conditional UI when the peptide correlation heatmap cannot be generated.
  output$pq_condition_no_peptide_heatmap <- shiny::renderUI({
    
    shiny::req(input$pq_body_condition_table_rows_selected, is.null(peptideHeatmapCondition()))
    
    shiny::tags$p(shiny::tags$h5(id = "pq_condition_no_peptide_heatmap_text", style = "display:inline", "Peptide correlation heatmap cannot be generated for this protein group."),
                  shiny::tags$h6(id = "pq_condition_no_peptide_heatmap_icon", style = "display:inline", icon("question-circle")),
                  shinyBS::bsTooltip(id = "pq_condition_no_peptide_heatmap_icon", title = "Too many missing values in peptide abundance vectors, which successively lead to correlation vectors with less than 2 finite values or a single peptide remained to that protein group.", placement = "right", trigger = "hover", options = list(container = "body")))
    
  })
  
  #Renders peptideHeatmapCondition() to condition tab.
  output$pq_condition_peptide_heatmap <- plotly::renderPlotly({
    
    shiny::req(peptideHeatmapCondition())
    
    return(peptideHeatmapCondition())
    
  })
  
  #Button to download peptide abundance correlation heatmap in static version.
  output$pq_condition_peptide_heatmap_download_cui <- shiny::renderUI({
    
    shiny::req(quantificationCondition(), probesAndWeights(), peptideHeatmapCondition())

    if(!is.null(peptideHeatmapCondition())){
      
      shiny::tags$p(shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_condition_peptide_heatmap_width", label = "Width: ", min = 100, max = 1000, value = 1000, step = 100, width = "250px")),
                    shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::numericInput("pq_condition_peptide_heatmap_height", label = "Height: ", min = 100, max = 1000, value = 500, step = 100, width = "200px")),
                    shiny::tags$div(style = "display: inline-block;margin: 25px 0px 0px 0px;", shiny::downloadButton("pq_condition_peptide_heatmap_download", "Download")),
                    shinyBS::bsTooltip(id = "pq_condition_peptide_heatmap_download", title = "Choose the width and height and download the above figure in PDF format.", placement = "right", trigger = "hover", options = list(container = "body")))
    }
    
  })
  
  #Save peptide abundance correlation heatmap in PDF format.
  output$pq_condition_peptide_heatmap_download <- shiny::downloadHandler(
    
    filename = function() {
      
      paste("PeptideAbundanceHeatmapCondition_", gsub("\\|","_",probesAndWeights()$protein.group[input$pq_body_condition_table_rows_selected]), "_",Sys.time(), ".pdf", collapse = "", sep = "")
      
    },
    content = function(file) {
      
      temp_name <- tempfile(pattern = "peptideHeatmap", fileext = ".html")
      
      p <- peptideHeatmapCondition()
      
      if(input$pq_condition_peptide_heatmap_width > 0 & input$pq_condition_peptide_heatmap_height > 0){
        
        p$width <- input$pq_condition_peptide_heatmap_width
        p$height <- input$pq_condition_peptide_heatmap_height
        
        htmlwidgets::saveWidget(p, temp_name)
        
        webshot2::webshot(url = temp_name, file = file)
        
        unlink(temp_name)
        unlink(paste(gsub(".html", "", temp_name), "_files", collapse = "", sep = ""), recursive = T)
        
        shiny::updateNumericInput(session, "pq_condition_peptide_heatmap_width", value = 1000)
        shiny::updateNumericInput(session, "pq_condition_peptide_heatmap_height", value = 500)
        
      }
      
    }
    
  )
  
  ########################### VISUALIZATION MODULES FOR CONDITIONS TAB END###########################
  
  output$pq_tabbox <- shiny::renderUI({
    
    shiny::req(input$pq_quantification_go, quantificationThreshold())
    
    if(input$cns_total_quantitative_columns != input$cns_total_conditions){
      
      shiny::fluidRow(
        shinydashboard::tabBox(title = "Protein Quantification",
                               id = "pq_tab", 
                               width = 12,
                               shiny::tabPanel("Per sample",
                                               shiny::tags$p(shiny::tags$h5(id = "pq__body_text1", style="display:inline;", shiny::tags$b("Protein quantification per sample.")),
                                                             shiny::tags$h6(id = "pq_body_q_icon1", style = "display:inline;", icon("question-circle"))),
                                               shinyBS::bsTooltip(id = "pq_body_q_icon1", title = "Change the peptide weight threshold by the slider below. The small plot on the left is the Global Correlation Index (GCI) of the protein expression profiles per peptide weight threshold value. Within trace (blue) corresponds to the mean correlation of samples within the same condition while between trace (orange) corresponds to the absolute mean correlation between different conditions. Smiley faces are used to indicate the peptide weight threshold that maximizes the within correlation.", placement = "right", trigger = "hover", options = list(container = "body")),
                                               shiny::tags$p(shiny::tags$div(style = "display: inline-block; vertical-align:top;", shiny::uiOutput("pq_body_weight_slider")),
                                                             shiny::tags$div(style = "display: inline-block; vertical-align:top;", plotly::plotlyOutput("pq_correlation_indicator", width = 280, height = 80)),
                                                             shiny::tags$div(style = "display: inline-block; vertical-align:top; margin: 0px 0px 0px 20px;", shiny::imageOutput("pq_smileyImage", height = 70, width = 70))),
                                               shiny::fluidRow(shinydashboard::box(id = "pq_box1", title = "Protein Quantification table", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary",
                                                                shinycssloaders::withSpinner(DT::dataTableOutput("pq_body_threshold_table"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                shiny::tags$style(type = "text/css", "#pq_threshold_table_download {border-color:#2A9C7F;}"),
                                                                shiny::column(width = 4, shiny::tags$div(style = "display: inline-block; vertical-align:top; margin-top:5px;", shiny::uiOutput("pq_threshold_table_download_button")), style = "padding-left: 0px;"),
                                                                shiny::column(width = 8, shiny::tags$div(style = "display: inline-block; vertical-align:top; margin-left:3px; margin-top:5px; float:right;", shiny::uiOutput("pq_total_sum_sample_download_button")),
                                                                                         shiny::tags$div(style = "display: inline-block; vertical-align:top; margin-top:5px; float:right;", shiny::uiOutput("pq_weighted_sum_sample_download_button")), style = "padding-right: 0px;"),
                                                                width = 12)),
                                               shiny::fluidRow(shinydashboard::box(id = "pq_box3", title = "Protein expression correlation heatmap", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary",
                                                                shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_clust_select", label = "Select clustesting method:", choices = c("Ward's", "Single linkage", "Complete linkage", "UPGMA", "WPGMA"), selected = "Complete linkage", multiple = F)),
                                                                shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_distance_select", label = "Select distance measure:", choices = c("Euclidean", "Maximum", "Manhattan", "Canberra", "Binary", "Minkowski"), selected = "Euclidean", multiple = F)),
                                                                shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_n_clusters", label = "Select number of clusters:", choices = c("Any", 1:input$cns_total_quantitative_columns), selected = "Any", multiple = F)),
                                                                shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_gradient", label = "Choose heatmap gradient:", choices = c("Gradient 1" = "1.png", "Gradient 2" = "2.png", "Gradient 3" = "3.png", "Gradient 4" = "4.png"), selected = "2.png", multiple = F, options = list(render = I("{ option: function(item, escape) { return '<div><img src=\"' + item.value + '\" width = 150 />' + escape(item.label) + '</div>' } }"))  )), 
                                                                shinycssloaders::withSpinner(plotly::plotlyOutput("pq_threshold_protein_heatmap"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                shiny::uiOutput("pq_threshold_protein_heatmap_download_cui"), width = 12)),
                                               shinyBS::bsCollapse(id = "pq_collapse1",
                                                                   shinyBS::bsCollapsePanel(title = "Protein group line plot", value = "pq_panel1", style = "primary",
                                                                    shinycssloaders::withSpinner(plotly::plotlyOutput("pq_threshold_lineplot"), type = 6, color = "#2A9C7F", size = 0.5) ,
                                                                    shiny::uiOutput("pq_threshold_lineplot_download_cui"))),
                                               shinyBS::bsCollapse(id = "pq_collapse3",
                                                                   shinyBS::bsCollapsePanel(title = "Protein group VIQoR plot", value = "pq_panel3", style = "primary",
                                                                    shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectInput("pq_threshold_select1", label = "Reference sample:", choices = NULL, selected = NULL)),
                                                                    shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectInput("pq_threshold_select2", label = "Testing sample:", choices = NULL, selected = NULL)),
                                                                    shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectInput("pq_threshold_enzyme", label = "Enzyme:", choices = c("None", "Trypsin", "Trypsin Strict", "Chymotrypsin (High)", "Chymotrypsin (Low)", "Pepsin (pH1.3)", "Pepsin (pH>2)", "Lys-C", "Arg-C"), selected = "None")),
                                                                    shinycssloaders::withSpinner(plotly::plotlyOutput("pq_threshold_VIQoR"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                    shiny::tags$style(".form-control.shiny-bound-input.shinyjs-resettable{border-radius:4px;}"),
                                                                    shiny::tags$style(".form-control:focus.shiny-bound-input.shinyjs-resettable{border-color:#66afe9;outline:0;box-shadow:0 0 8px rgb(102,175,233,0.6)}"),
                                                                    shiny::uiOutput("pq_threshold_VIQoR_download_cui"))),
                                               shinyBS::bsCollapse(id = "pq_collapse5",
                                                                   shinyBS::bsCollapsePanel(title = "Peptide abundance correlation heatmap", value = "pq_panel5", style = "primary",
                                                                    shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_clust_select_pep", label = "Select clustesting method:", choices = c("Ward's", "Single linkage", "Complete linkage", "UPGMA", "WPGMA"), selected = "Complete linkage", multiple = F)),
                                                                    shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_distance_select_pep", label = "Select distance measure:", choices = c("Euclidean", "Maximum", "Manhattan", "Canberra", "Binary", "Minkowski"), selected = "Euclidean", multiple = F)),
                                                                    shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_n_clusters_pep", label = "Select number of clusters:", choices = c("Any", 1), selected = "Any", multiple = F)),
                                                                    shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_gradient_pep", label = "Choose heatmap gradient:", choices = c("Gradient 1" = "1.png", "Gradient 2" = "2.png", "Gradient 3" = "3.png", "Gradient 4" = "4.png"), selected = "2.png", multiple = F, options = list(render = I("{ option: function(item, escape) { return '<div><img src=\"' + item.value + '\" width = 150 />' + escape(item.label) + '</div>' } }"))  )), 
                                                                    shiny::uiOutput("pq_threshold_no_peptide_heatmap"),
                                                                    shinycssloaders::withSpinner(plotly::plotlyOutput("pq_threshold_peptide_heatmap"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                    shiny::uiOutput("pq_threshold_peptide_heatmap_download_cui")))
                                              ),
                               if(input$cns_total_conditions != 1){
                                  shiny::tabPanel("Per Condition", 
                                                  shiny::tags$p(shiny::tags$h5(id = "pq__body_text2", style="display:inline;", shiny::tags$b("Protein quantification per condition.")),
                                                                shiny::tags$h6(id = "pq_body_q_icon2", style = "display:inline;", icon("question-circle"))),
                                                  shinyBS::bsTooltip(id = "pq_body_q_icon2", title = "Protein expression is summarized per condition by average, for the selected peptide weight threshold.", placement = "right", trigger = "hover", options = list(container = "body")),
                                                  shiny::fluidRow(shinydashboard::box(id = "pq_box2", title = "Protein Quantification table", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary",
                                                                  shinycssloaders::withSpinner(DT::dataTableOutput("pq_body_condition_table"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                  shiny::tags$style(type = "text/css", "#pq_condition_table_download {border-color:#2A9C7F;}"),
                                                                  shiny::column(width = 4, shiny::tags$div(style = "display: inline-block; vertical-align:top; margin-top:5px;", shiny::uiOutput("pq_condition_table_download_button")), style = "padding-left: 0px;"), 
                                                                  shiny::column(width = 8, shiny::tags$div(style = "display: inline-block; vertical-align:top; margin-left:3px; margin-top:5px; float:right;", shiny::uiOutput("pq_weighted_sum_condition_download_button")),
                                                                                           shiny::tags$div(style = "display: inline-block; vertical-align:top; margin-top:5px; float:right;", shiny::uiOutput("pq_total_sum_condition_download_button")), style = "padding-right: 0px;"),
                                                                  width = 12)),
                                                  shiny::fluidRow(shinydashboard::box(id = "pq_box4", title = "Protein expression correlation heatmap", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary",
                                                                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_condition_clust_select", label = "Select clustesting method:", choices = c("Ward's", "Single linkage", "Complete linkage", "UPGMA", "WPGMA"), selected = "Complete linkage", multiple = F)),
                                                                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_condition_distance_select", label = "Select distance measure:", choices = c("Euclidean", "Maximum", "Manhattan", "Canberra","Binary", "Minkowski"), selected = "Euclidean", multiple = F)),
                                                                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_condition_n_clusters", label = "Select number of clusters:", choices = c("Any", 1, 2), selected = 2, multiple = F)),
                                                                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_condition_gradient", label = "Choose heatmap gradient:", choices = c("Gradient 1" = "1.png", "Gradient 2" = "2.png", "Gradient 3" = "3.png", "Gradient 4" = "4.png"), selected = "2.png", multiple = F, options = list(render = I("{ option: function(item, escape) { return '<div><img src=\"' + item.value + '\" width = 150 />' + escape(item.label) + '</div>' } }"))  )), 
                                                                  shinycssloaders::withSpinner(plotly::plotlyOutput("pq_condition_protein_heatmap"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                  shiny::uiOutput("pq_condition_protein_heatmap_download_cui"), width = 12)),
                                                  shinyBS::bsCollapse(id = "pq_collapse2",
                                                                      shinyBS::bsCollapsePanel(title = "Protein group line plot", value = "pq_panel2", style = "primary",
                                                                      shinycssloaders::withSpinner(plotly::plotlyOutput("pq_condition_lineplot"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                      shiny::uiOutput("pq_condition_lineplot_download_cui"))),
                                                  shinyBS::bsCollapse(id = "pq_collapse4",
                                                                      shinyBS::bsCollapsePanel(title = "Protein group VIQoR plot", value = "pq_panel4", style = "primary",
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectInput("pq_condition_select1", label = "Reference condition:", choices = NULL, selected = NULL)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectInput("pq_condition_select2", label = "Testing condition:", choices = NULL, selected = NULL)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectInput("pq_condition_enzyme", label = "Enzyme:", choices = c("None", "Trypsin", "Trypsin Strict", "Chymotrypsin (High)", "Chymotrypsin (Low)", "Pepsin (pH1.3)", "Pepsin (pH>2)", "Lys-C", "Arg-C"), selected = "None")),
                                                                      shinycssloaders::withSpinner(plotly::plotlyOutput("pq_condition_VIQoR"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                      shiny::uiOutput("pq_condition_VIQoR_download_cui"))),
                                                  shinyBS::bsCollapse(id = "pq_collapse6",
                                                                      shinyBS::bsCollapsePanel(title = "Peptide abundance correlation heatmap", value = "pq_panel6", style = "primary",
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_condition_clust_select_pep", label = "Select clustesting method:", choices = c("Ward's", "Single linkage", "Complete linkage", "UPGMA", "WPGMA"), selected = "Complete linkage", multiple = F)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_condition_distance_select_pep", label = "Select distance measure:", choices = c("Euclidean", "Maximum", "Manhattan", "Canberra","Binary", "Minkowski"), selected = "Euclidean", multiple = F)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_condition_n_clusters_pep", label = "Select number of clusters:", choices = c("Any", 1, 2), selected = 2, multiple = F)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_condition_gradient_pep",label = "Choose heatmap gradient:", choices = c("Gradient 1" = "1.png", "Gradient 2" = "2.png", "Gradient 3" = "3.png", "Gradient 4" = "4.png"), selected = "2.png", multiple = F, options = list(render = I("{ option: function(item, escape) { return '<div><img src=\"' + item.value + '\" width = 150 />' + escape(item.label) + '</div>' } }"))  )), 
                                                                      shiny::uiOutput("pq_condition_no_peptide_heatmap"),
                                                                      shinycssloaders::withSpinner(plotly::plotlyOutput("pq_condition_peptide_heatmap"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                      shiny::uiOutput("pq_condition_peptide_heatmap_download_cui")))
                                                )
                                              }
                              )
                      )
      
    } else {
      
        shiny::fluidRow(
          shinydashboard::tabBox(title = "Protein Quantification",
                                 id = "pq_tab", 
                                 width = 12,
                                 shiny::tabPanel("Per condition",
                                                 shiny::tags$p(shiny::tags$h5(id = "pq__body_text1", style="display:inline;", shiny::tags$b("Protein quantification per condition.")),
                                                               shiny::tags$h6(id = "pq_body_q_icon1", style = "display:inline;", icon("question-circle"))),
                                                 shinyBS::bsTooltip(id = "pq_body_q_icon1", title = "Change the peptide weight threshold by the slider below. The small plot on the left is the Global Correlation Index (GCI) of the protein expression profiles per peptide weight threshold value. The between trace (blue) corresponds to the absolute mean correlation between different conditions.", placement = "right", trigger = "hover", options = list(container = "body")),
                                                 shiny::tags$p(shiny::tags$div(style = "display: inline-block; vertical-align:top;", shiny::uiOutput("pq_body_weight_slider")),
                                                               shiny::tags$div(style = "display: inline-block; vertical-align:top;", plotly::plotlyOutput("pq_correlation_indicator", width = 280, height = 80)),
                                                               shiny::tags$div(style = "display: inline-block; vertical-align:top; margin: 0px 0px 0px 20px;", shiny::imageOutput("pq_smileyImage", height = 70, width = 70))),
                                                 shiny::fluidRow(shinydashboard::box(id = "pq_box1", title = "Protein Quantification table", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary",
                                                                 shinycssloaders::withSpinner(DT::dataTableOutput("pq_body_threshold_table"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                 shiny::tags$style(type = "text/css", "#pq_threshold_table_download {border-color:#2A9C7F;}"),
                                                                 shiny::column(width = 4, shiny::tags$div(style = "display: inline-block; vertical-align:top; margin-top:5px;", shiny::uiOutput("pq_threshold_table_download_button")), style = "padding-left: 0px;"),
                                                                 shiny::column(width = 8, shiny::tags$div(style = "display: inline-block; vertical-align:top; margin-left:3px; margin-top:5px; float:right;", shiny::uiOutput("pq_total_sum_sample_download_button")),
                                                                               shiny::tags$div(style = "display: inline-block; vertical-align:top; margin-top:5px; float:right;", shiny::uiOutput("pq_weighted_sum_sample_download_button")), style = "padding-right: 0px;"),
                                                                 width = 12)),
                                                 shiny::fluidRow(shinydashboard::box(id = "pq_box3", title = "Protein expression correlation heatmap", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE, status = "primary",
                                                                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_clust_select", label = "Select clustesting method:", choices = c("Ward's", "Single linkage", "Complete linkage", "UPGMA", "WPGMA"), selected = "Complete linkage", multiple = F)),
                                                                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_distance_select", label = "Select distance measure:", choices = c("Euclidean", "Maximum", "Manhattan", "Canberra", "Binary", "Minkowski"), selected = "Euclidean", multiple = F)),
                                                                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_n_clusters", label = "Select number of clusters:", choices = c("Any", 1), selected = "Any", multiple = F)),
                                                                  shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_gradient", label = "Choose heatmap gradient:", choices = c("Gradient 1" = "1.png", "Gradient 2" = "2.png", "Gradient 3" = "3.png", "Gradient 4" = "4.png"), selected = "2.png", multiple = F, options = list(render = I("{ option: function(item, escape) { return '<div><img src=\"' + item.value + '\" width = 150 />' + escape(item.label) + '</div>' } }"))  )), 
                                                                  shinycssloaders::withSpinner(plotly::plotlyOutput("pq_threshold_protein_heatmap"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                  shiny::uiOutput("pq_threshold_protein_heatmap_download_cui"), width = 12)),
                                                 shinyBS::bsCollapse(id = "pq_collapse1",
                                                                     shinyBS::bsCollapsePanel(title = "Protein group line plot", value = "pq_panel1", style = "primary",
                                                                      shinycssloaders::withSpinner(plotly::plotlyOutput("pq_threshold_lineplot"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                      shiny::uiOutput("pq_threshold_lineplot_download_cui"))),
                                                 shinyBS::bsCollapse(id = "pq_collapse3",
                                                                     shinyBS::bsCollapsePanel(title = "Protein group VIQoR plot", value = "pq_panel3", style = "primary",
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectInput("pq_threshold_select1", label = "Reference sample:", choices = NULL, selected = NULL)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectInput("pq_threshold_select2", label = "Testing sample:", choices = NULL, selected = NULL)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectInput("pq_threshold_enzyme", label = "Enzyme:", choices = c("None", "Trypsin", "Trypsin Strict", "Chymotrypsin (High)", "Chymotrypsin (Low)", "Pepsin (pH1.3)", "Pepsin (pH>2)", "Lys-C", "Arg-C"), selected = "None")),
                                                                      shinycssloaders::withSpinner(plotly::plotlyOutput("pq_threshold_VIQoR"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                      shiny::tags$style(".form-control.shiny-bound-input.shinyjs-resettable{border-radius:4px;}"),
                                                                      shiny::tags$style(".form-control:focus.shiny-bound-input.shinyjs-resettable{border-color:#66afe9;outline:0;box-shadow:0 0 8px rgb(102,175,233,0.6)}"),
                                                                      shiny::uiOutput("pq_threshold_VIQoR_download_cui"))),
                                                 shinyBS::bsCollapse(id = "pq_collapse5",
                                                                     shinyBS::bsCollapsePanel(title = "Peptide abundance correlation heatmap", value = "pq_panel5", style = "primary",
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_clust_select_pep", label = "Select clustesting method:", choices = c("Ward's", "Single linkage", "Complete linkage", "UPGMA", "WPGMA"), selected = "Complete linkage", multiple = F)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_distance_select_pep", label = "Select distance measure:", choices = c("Euclidean", "Maximum", "Manhattan", "Canberra", "Binary", "Minkowski"), selected = "Euclidean", multiple = F)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_n_clusters_pep", label = "Select number of clusters:", choices = c("Any", 1), selected = "Any", multiple = F)),
                                                                      shiny::tags$div(style = "display: inline-block;vertical-align:top; width: 200px;", shiny::selectizeInput("pq_threshold_gradient_pep", label = "Choose heatmap gradient:", choices = c("Gradient 1" = "1.png", "Gradient 2" = "2.png", "Gradient 3" = "3.png", "Gradient 4" = "4.png"), selected = "2.png", multiple = F, options = list(render = I("{ option: function(item, escape) { return '<div><img src=\"' + item.value + '\" width = 150 />' + escape(item.label) + '</div>' } }"))  )), 
                                                                      shiny::uiOutput("pq_threshold_no_peptide_heatmap"),
                                                                      shinycssloaders::withSpinner(plotly::plotlyOutput("pq_threshold_peptide_heatmap"), type = 6, color = "#2A9C7F", size = 0.5),
                                                                      shiny::uiOutput("pq_threshold_peptide_heatmap_download_cui")))
                                                 )
                                )
                        )
    }
    
  })
  
  ########################### PROTEIN QUANTIFICATION TAB END ###########################
  
  palette <- c('#0911e8', '#f50a0a', '#822205', '#8f6606', '#839c06', '#51cf08', '#09e854', '#069c5b', '#08c8cf', '#056182', 
              '#062a8f', '#3b0582', '#a407b5', '#cf0886', '#f50a31', '#db0909', '#f5770a', '#e8d909', '#6d8205', '#42a807', 
              '#08cf4a', '#05824c', '#07afb5', '#0979e8', '#9e09e8', '#8d069c', '#b50775', '#db092c', '#c20808', '#db6b09', 
              '#c2b508', '#97f50a', '#388f06', '#07b541', '#0af5c6', '#06979c', '#0865c2', '#080fcf', '#8c08cf', '#760582', 
              '#9c0665', '#c20827', '#a80707', '#b55807', '#a89e07', '#87db09', '#21c208', '#069c38', '#09dbb1', '#057e82', 
              '#0758a8', '#060a8f', '#7b07b5', '#e809ca', '#820554', '#9c061f', '#8f0606', '#8f4606', '#8f8606', '#77c208', 
              '#1a9c06', '#05822f', '#08c29d', '#0ab6f5', '#064a8f', '#2707a8', '#6a069c', '#c208a9', '#e80962', '#82051a', 
              '#e83d09', '#e8a509', '#c3e809', '#68a807', '#0af521', '#09e888', '#07a888', '#089acf', '#0a48f5', '#6a09e8', 
              '#580582', '#9c0688', '#c20852', '#c23308', '#c28a08', '#aecf08', '#588f06', '#09db1e', '#08cf79', '#068f73', 
              '#0787b5', '#0941db', '#5808c2', '#d209e8', '#820571', '#9c0642', '#a82c07', '#a87807', '#98b507', '#5be809', 
              '#068f13', '#07b56a', '#09e1e8', '#06749c', '#0736b5', '#47069c', '#bb08cf', '#e80996', '#820537')
  
  modifications.palette <- c("#f58231", "#911eb4", "#e6194B", "#ffe119", "#4363d8", "#3cb44b", "#bfef45", "#800000", "#aaffc3", "#9A6324")
  
}


shinyApp(ui, server)





