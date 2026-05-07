library(shiny)
library(survival)
library(bslib)


bridge_list <- readRDS("bridge_list.rds")

ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#2C3E50",
    base_font = font_google("Source Sans 3")
  ),
  
  div(
    style = "max-width: 1050px; margin: 0 auto; padding: 30px 20px;",
    
    div(
      style = "text-align: center; margin-bottom: 30px;",
      h2("BRIDGE Prognostic Index for DLBCL"),
      p(
        "Evaluate prognostic risk group and predicted 5-year overall survival probability ",
        "for patients with diffuse large B-cell lymphoma.",
        style = "font-size: 17px; color: #555;"
      )
    ),
    
    card(
      card_header(tags$b("Patient Characteristics")),
      
      layout_columns(
        col_widths = c(4, 4, 4),
        
        selectInput(
          inputId = "age_cat",
          label = "Age",
          choices = c(
            "Select age group" = "",
            "<=40" = "<=40",
            "41-60" = "41-60",
            "61-75" = "61-75",
            ">75" = ">75"
          ),
          selected = ""
        ),
        
        selectInput(
          inputId = "stage",
          label = "Ann Arbor stage",
          choices = c(
            "Select Ann Arbor stage" = "",
            "Stage I/II" = "I/II",
            "Stage III/IV" = "III/IV"
          ),
          selected = ""
        ),
        
        selectInput(
          inputId = "ldh",
          label = "Serum lactate dehydrogenase (LDH)",
          choices = c(
            "Select LDH status" = "",
            "<= upper limit of normal" = "<= upper limit of normal",
            "> upper limit of normal" = "> upper limit of normal"
          ),
          selected = ""
        )
      ),
      
      layout_columns(
        col_widths = c(6, 6),
        
        selectInput(
          inputId = "ecog",
          label = "ECOG performance status",
          choices = c(
            "Select ECOG PS" = "",
            "0-1" = "0-1",
            "2-4" = "2-4"
          ),
          selected = ""
        ),
        
        selectInput(
          inputId = "extranodal",
          label = "Number of extranodal sites",
          choices = c(
            "Select extranodal site category" = "",
            "0 or 1" = "0 or 1",
            "> 1" = "> 1"
          ),
          selected = ""
        )
      ),
      
      div(
        style = "text-align: center; margin-top: 20px;",
        actionButton(
          inputId = "run_prediction",
          label = "Run prediction",
          class = "btn-primary btn-lg"
        )
      )
    ),
    
    br(),
    
    card(
      card_header(tags$b("Prediction Result")),
      uiOutput("risk_result")
    )
  )
)

server <- function(input, output, session) {
  
  make_new_patient <- function() {
    
    dat <- data.frame(
      age_cat = input$age_cat,
      stage = input$stage,
      LDH = input$ldh,
      ECOG = input$ecog,
      extranodal = input$extranodal,
      stringsAsFactors = FALSE
    )
    
    ## Age dummy variables
    ## Reference group: 61-75
    dat$age_cat_le40 <- ifelse(dat$age_cat == "<=40", 1, 0)
    dat$age_cat_41_60 <- ifelse(dat$age_cat == "41-60", 1, 0)
    dat$age_cat_gt75 <- ifelse(dat$age_cat == ">75", 1, 0)
    
    ## Other IPI variables
    dat$stage_IPI <- ifelse(dat$stage == "III/IV", 1, 0)
    dat$LDH.IPI <- ifelse(dat$LDH == "> upper limit of normal", 1, 0)
    dat$ECOG.IPI <- ifelse(dat$ECOG == "2-4", 1, 0)
    dat$multiple.extranodal.IPI <- ifelse(dat$extranodal == "> 1", 1, 0)
    
    dat
  }
  
  get_beta_vector <- function() {
    
    beta_use <- bridge_list$beta
    
    ## Convert matrix/data.frame beta to vector if needed
    if (is.matrix(beta_use) || is.data.frame(beta_use)) {
      beta_names <- rownames(beta_use)
      
      if (is.null(beta_names)) {
        beta_names <- colnames(beta_use)
      }
      
      beta_use <- as.numeric(beta_use)
      
      if (!is.null(beta_names) && length(beta_names) == length(beta_use)) {
        names(beta_use) <- beta_names
      }
    }
    
    ## If beta has no names, assign X_var names by position
    if (is.null(names(beta_use))) {
      if (length(beta_use) != length(bridge_list$X_var)) {
        stop("beta has no names and its length does not match bridge_list$X_var.")
      }
      
      names(beta_use) <- bridge_list$X_var
    }
    
    if (!all(bridge_list$X_var %in% names(beta_use))) {
      stop("The names of beta must include all variables in bridge_list$X_var.")
    }
    
    beta_use
  }
  
  prediction <- eventReactive(input$run_prediction, {
    
    required_inputs <- c(
      input$age_cat,
      input$stage,
      input$ldh,
      input$ecog,
      input$extranodal
    )
    
    validate(
      need(
        all(required_inputs != ""),
        "Please complete all five patient characteristics before running the prediction."
      )
    )
    
    dat <- make_new_patient()
    
    if (!all(bridge_list$X_var %in% colnames(dat))) {
      stop("Some X_var variables are missing from the new patient data.")
    }
    
    beta_use <- get_beta_vector()
    
    ## Compute offset from X variables
    dat$offset <- as.numeric(
      as.matrix(dat[, bridge_list$X_var, drop = FALSE]) %*%
        beta_use[bridge_list$X_var]
    )
    
    ## Newdata for Cox model: Z variables plus offset
    pred_dat <- dat[, c(bridge_list$Z_var, "offset"), drop = FALSE]
    
    sf <- survfit(
      bridge_list$reddy.cox,
      newdata = pred_dat
    )
    
    surv_5yr <- summary(sf, times = 5, extend = TRUE)$surv
    surv_5yr <- as.numeric(surv_5yr)
    
    risk_group <- ifelse(
      surv_5yr <= bridge_list$cutoff[1],
      "High",
      ifelse(
        surv_5yr <= bridge_list$cutoff[2],
        "Intermediate",
        "Low"
      )
    )
    
    list(
      surv_5yr = surv_5yr,
      risk_group = risk_group
    )
  })
  
  output$risk_result <- renderUI({
    
    if (input$run_prediction == 0) {
      return(
        div(
          style = "color: #666; font-size: 16px;",
          "Select all patient characteristics and click ",
          tags$b("Run prediction"),
          " to calculate the BRIDGE Prognostic Index risk group."
        )
      )
    }
    
    pred <- prediction()
    
    group_color <- switch(
      pred$risk_group,
      "Low" = "#198754",
      "Intermediate" = "#F39C12",
      "High" = "#C0392B"
    )
    
    div(
      style = "padding: 15px 5px;",
      
      div(
        style = paste0(
          "font-size: 30px; font-weight: 700; color: ",
          group_color,
          "; margin-bottom: 12px;"
        ),
        pred$risk_group,
        " Risk"
      ),
      
      div(
        style = "font-size: 21px; margin-bottom: 8px;",
        "Predicted 5-year overall survival probability: ",
        tags$b(paste0(round(pred$surv_5yr * 100, 1), "%"))
      ),
      
      div(
        style = "font-size: 15px; color: #666;",
        "Risk group was assigned using the BRIDGE Prognostic Index cutoffs at 64.9% (40th percentile) and 37.5% (10th percentile) defined the low (5-year OS > 64.9%), intermediate (37.5%-64.9%), and high (≤ 37.5%) BPI groups."
      )
    )
  })
}

shinyApp(ui = ui, server = server)