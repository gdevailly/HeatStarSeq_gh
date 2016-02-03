library(shiny)
library(shinyjs)
library(plotly)
library(magrittr)

# loading data painfully generated before
load("data/encode_rnaseq_preload.RData")
load("data/bgee_human_preload.RData")
load("data/encode_mouse_rnaseq_preload.RData")
load("data/bgee_mouse_preload.RData")
