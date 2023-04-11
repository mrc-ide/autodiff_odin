#Convert between the step of the adjoint model (reverse) and the main model (forward)
main_step <- total_steps - step - 1

#pull the states of the main model
time <- main_states[1,main_step+1]
S <- main_states[2,main_step+1]
R <- main_states[3,main_step+1]
I <- main_states[4,main_step+1]
cases_cumul <- main_states[5,main_step+1]
cases_inc <- main_states[6,main_step+1]

#Recalculate the intermediate variables between the state[main_step] and state[main_step+1]
p_IR <- 1 - exp(-(gamma))
freq <- user(4)
dt <- 1.0 / freq

N <- S + I + R
p_inf <- beta * I / N
p_SI <- 1 - exp(-(p_inf))
n_SI <- S * p_SI * dt
n_IR <- I * p_IR * dt

#Now starts the back propagation throughout the code ladder
#This is the "smart" bit to be built
#This should be relatively easily done
#Here I put some code that I have used to "manually" derive the equations
#First it builds the graph of dependencies in the main odin model
#
#parsed_model <- jsonlite::fromJSON(odin::odin_parse("models/sir_4_AD.R"))
# construct_param_tree <- function(parsed_model){
#   parameter_graph <- NULL
#   for(n in seq_along(parsed_model$equations$name)){
#     if(!is.null(parsed_model$equations$depends$variables[n][[1]])){
#       for(p in eval(parsed_model$equations$depends$variables[n][[1]]))
#         parameter_graph <- rbind(parameter_graph,
#                                  c(p,parsed_model$equations$name[n]))
#     }
#   }
#   igraph::graph_from_edgelist(parameter_graph)
# }
#param_graph <- construct_param_tree(parsed_model)
#names(param_graph[1,])[which(param_graph["n_IR",]==1)]
##>this returns "update_I" "update_R"
##>then we can recover from the main odin model the equations for update(I) and update(R)
##>And do symbolic differentiaton (believe base R can do it)
##>the RHS of the update(I) eq is "I + n_SI - n_IR"
##>We can run
#ff <- parse(text="I + n_SI - n_IR")
#D(ff, "n_IR")
##>This returns -1
##>We have the first element of the back-propagation for adj_n_IR:
##> (-1)*adj_I which combines the partial derivative with the child value
##> this works in reverse than the main odin model, so back-propagate the change in the Log-likelihood
##> We can do the same for the other "child" of n_IR
##> the RHS of the update(R) is "R + n_IR"
##> the derivation returns 1
##> so the adjoint equation of adj_n_IR is

adj_n_IR <- -adj_I + adj_R

#the same process can be applied for the other intermediate
#equations, giving:

adj_n_SI <- adj_cases_cumul + adj_cases_inc + adj_I - adj_S
adj_p_SI <- S * dt * adj_n_SI
adj_p_inf <- exp(-(p_inf)) * adj_p_SI
adj_N <- -(beta * I/N^2) * adj_p_inf

adj_p_IR <- I * dt * adj_n_IR

#adjoint variables: propagates feedback from
update(adj_time) <- 0
update(adj_S) <- adj_N + p_SI * dt * adj_n_SI + adj_S
update(adj_R) <- adj_N + adj_R
update(adj_I) <- adj_N + p_IR * dt * adj_n_IR + beta/N * adj_p_inf + adj_I
update(adj_cases_cumul) <- adj_cases_cumul
update(adj_cases_inc) <- if ((main_step!=0)&&((main_step) %% freq == 0)) data_input[main_step]/cases_inc - 1 else adj_cases_inc

#adjoint_parameters: accumulates feedback
update(adj_beta) <- adj_beta + I / N * adj_p_inf
update(adj_gamma) <- adj_gamma + exp(-(gamma)) * adj_p_IR
update(adj_I0) <- if(main_step == 0) adj_N + p_IR * dt * adj_n_IR + beta/N * adj_p_inf + adj_I else 0

initial(adj_time) <- 0
initial(adj_S) <- 0
initial(adj_R) <- 0
initial(adj_I) <- 0
initial(adj_cases_cumul) <- 0
initial(adj_cases_inc) <- if((total_steps) %% freq == 0) data_input[total_steps]/main_states[6,total_steps+1] - 1 else 0

main_states[,] <- user()
dim(main_states) <- user()

data_input[] <- user()
dim(data_input) <- user()

total_steps <- user()

initial(adj_beta) <- 0
initial(adj_gamma) <- 0
initial(adj_I0) <- 0

beta <- user(0.2)
gamma <- user(0.1)
I0 <- user(10)
