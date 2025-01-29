########## R CODE FOR SHINY APP FOR REPORTING TIME SIMULATIONS #################
### Packages required ###
pkgs <- c(
  "survival",
  "tinytex",
  "tidyverse",
  "dplyr",
  "flexsurv",
  "simsurv",
  "shiny",
  "plotly",
  "ggplot2",
  "bslib"
)

# Check if the package is installed, if yes load, if no: install + load 
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

rm(pkgs)

source( "code/compile.R")


############################## R SHINY USER INTERFACE ##########################
ui <- page_fluid(
  
  # Application title
  titlePanel("Early Report Sample Size Simulation"),
    
  # CARDS FOR DIFFERENT TOPICS:
  layout_columns(
    card(
      card_header("Trial Design"),
      numericInput("t_enroll", "Enrollment Period (years):", value = 1, min = 0.1),
      numericInput("t_outcome", "Outcome Time (years):", value = 1, min = 0.1),
      numericInput("n_patients", "Target Sample Size:", value = 100, min = 1, step = 1),
      numericInput("zero_enroll", "Patients enrolled at time 0:", value = 0, min = 0),
      checkboxInput("treat", "Include Treatment Effect?", value = FALSE),
      # Conditional panel to display the hazard ratio input only if treatment is selected
      conditionalPanel(
        condition = "input.treat == true",
        numericInput("hr", "Hazard Ratio:", value = 1, min = 0.01),
        numericInput("pat2", "Patients in treatment group", value = 100, min = 10)
      ),
    ),
    
    card(
      card_header("Censoring Mechanism (WIP)"),
      checkboxInput("censor", "Include Censoring Hazard?", value = FALSE),
      conditionalPanel(
        condition = "input.censor == true",
        radioButtons("dist.c", "Distribution:", 
                     choices = c("Weibull" = "weibull", "Exponential" = "exponential", "Gompertz" = "gompertz")),
        numericInput("scale.c", "Scale:", value = 1, min = 0.1 ), 
        conditionalPanel(
          condition = "input.dist.c != 'exponential' ",
          numericInput("shape.c", "Shape:", value = 1, min = 0.1)
        ),
      ),
    ),
  
    card(
      card_header("Cumulative Baseline Hazard"),
      layout_sidebar(
        sidebar = sidebar(
          radioButtons("dist", "Distribution:", 
                      choices = c("Weibull" = "weibull", "Exponential" = "exponential", "Gompertz" = "gompertz")),
          numericInput("scale", "Scale:", value = 1, min = 0.1 ), 
          conditionalPanel(
            condition = "input.dist != 'exponential' ",
            numericInput("shape", "Shape:", value = 1, min = 0.1)
          ),
        ),
        plotOutput("cumHaz")
      )
    ),
    
    card(
      card_header( "Simulation and plot settings" ),
      actionButton("run_simulation", "Run Simulation"),
      numericInput("iter", "Iterations:", value = 10, min = 1),
      numericInput("points", "Points:", value = 20, min = 2),
      textInput("xlab", "X-axis Label:", value = "Reporting time after start enrollment"),
      textInput("ylab", "Y-axis Label:", value = "Percentage of total sample size"),
      textInput("title", "Plot Title:", value = "Effective sample size as percentage over time"),
      checkboxInput("plotly", "Interactive plot (WIP)", value = TRUE),
      checkboxInput("mod", "Modified Effective N", value = TRUE),
      checkboxInput("conf", "Show confidence interval", value = FALSE),
      conditionalPanel(
        condition = "input.conf == TRUE",
        sliderInput("conf.int", "Percentage Confidence", value = 95, min = 0, max = 100)
      )
    ),
    
    # Main panel for displaying the plot
    card(
      card_header("Simulation results"),
      plotOutput("simPlot")
    ),
    
    # Single instance
    card(
      card_header( "Results at time t (WIP)")
    ),
    
    card(
      card_header( "Here will be output" )
    ),
    
    col_widths = c(2, 2, 7, 3, 8, 3, 8)
  )
)



############################## R SHINY SERVER BACKEND ##########################
# Define the server logic to run the simulation
server <- function(input, output) {
  
  # Reactive expression to generate the plot based on user input
  observeEvent(input$run_simulation, {
    output$simPlot <- renderPlot({
      # Generate X depending on treatment
      if( input$treat ){
        x = data.frame( id = 1:(input$n_patients+input$pat2), 
                        "treat" = c(rep(0, input$n_patients), rep(1, input$pat2)))
        betas <- c( "treat" = log(input$hr))
      }
      else{ 
        x = data.frame(id = 1:input$n_patients)
        betas = NULL
      }
      
      # Call the early_report_ss function with user inputs
      plot <- trial_report(dist = input$dist, 
                      gammas = input$shape,
                      lambdas = input$scale^(-input$shape),
                      x = x,
                      zero_enroll = input$zero_enroll, 
                      t_enroll = input$t_enroll, 
                      t_outcome = input$t_outcome, 
                      iter = input$iter, 
                      points = input$points, 
                      xlab = input$xlab, 
                      ylab = input$ylab, 
                      title = input$title,
                      betas = betas,
                      mod = input$mod,
                      conf.int = input$conf.int,
                      conf = input$conf
      )
      return(plot)
    })
  })
  
  output$cumHaz <- renderPlot({
    # Sim the data
    custom_colors <- c("Baseline" = "#C2666B", 
                       "Treatment" = "#c6aa2c")
    
    if( input$dist == "weibull" ) data = -log( 1- pweibull(seq(0, input$t_outcome, length.out = 101), shape = input$shape, scale = input$scale) )
    else if( input$dist == "exponential") data = -log( 1- pweibull(seq(0, input$t_outcome, length.out = 101), shape = 1, scale = input$scale) )
    else if( input$dist == "gompertz" ) data = -log( 1- pgompertz(seq(0, input$t_outcome, length.out = 101), shape = input$shape, rate = input$scale) )

    
    df <- data.frame(
      x = seq(0, input$t_outcome, length.out = 101),
      y = data,
      group = "Baseline"
      )
    if( input$treat ) df <- rbind( df, 
                                   data.frame(
              x = seq(0, input$t_outcome, length.out = 101),
              y = input$hr*data,
              group = "Treatment"
            ))
            
    ggplot(df, aes(x = x, y = y, color = group)) +
      geom_line( linewidth = 1.2 ) +
      xlim(0, input$t_outcome) +
      ylim(0, 2) +
      labs(title = "Cumulative Hazard",
           x = "Time (years)",
           y = "Cumulative hazard") +
      theme_minimal()
  })
}

############################## RUN R SHINY ##########################
shinyApp(ui = ui, server = server)
