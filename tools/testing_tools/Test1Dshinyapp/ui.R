library(dygraphs)
library(shiny)

shinyUI(

fluidPage(
		
		titlePanel("GEOtop Test Cases (1D)"),
		
		column(4, wellPanel(
						uiOutput("sim_names_ui"),
						uiOutput("sim_kws_ui"),
						uiOutput("sim_layer_ui")
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