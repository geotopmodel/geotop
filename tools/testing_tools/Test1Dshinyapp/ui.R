library(dygraphs)
library(shiny)
library(shinyjs)

shinyUI(

fluidPage(
	
#		div(
#				id = "loading-content",
#				h2("Loading...")
#		),
		titlePanel("GEOtop Test Cases (1D)"),
		
		column(4, wellPanel(
						uiOutput("sim_names_ui"),
						uiOutput("sim_kws_ui"),
						uiOutput("sim_layer_ui"),
						checkboxInput(inputId = "latest",label = strong("latest"),value = TRUE),
						checkboxInput(inputId = "SE27XX",label = strong("SE27XX"),value = TRUE),
						checkboxInput(inputId = "METEOIO_OFF",label = strong("METEOIO-OFF"),value = TRUE),
						checkboxInput(inputId = "METEOIO_ON",label = strong("METEOIO-ON"),value = FALSE),
						sliderInput(inputId="digitsAfterDecimal", label=strong("Digits After Decimal:"),min=0, max=20, value=5)	
								
#uiOutput("sim_vars_ui")
						#selectInput('var', 'Variable',uiOutput("sim_vars"))
						#sliderInput("n", "Number of points:",
						#		min = 10, max = 200, value = 50, step = 10)
		#		)
				)),
		
		column(5,
				#"The plot below will be not displayed when the slider value",
				#"is less than 50.",
				"Values:",
				"",
				"",
				# With the conditionalPanel, the condition is a JavaScript
				# expression. In these expressions, input values like
				# input$n are accessed with dots, as in input.n
				#showPanel("input.n >= 50",
						dygraphOutput("dygraph"), ##, height = 300)
						"",
						"",
				 "residuals from SE27XX:",
				 "",
				 "",
				        dygraphOutput("dygraph_res")
				#)
		)

))