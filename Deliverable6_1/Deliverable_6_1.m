clc
close all
clear all
import casadi.*

quad = Quad();
CTRL= ctrl_NMPC(quad);

sim = quad.sim(CTRL);
quad.plot(sim);