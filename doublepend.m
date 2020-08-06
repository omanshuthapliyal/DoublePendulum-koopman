function [T, Y] = doublepend(tspan, initCond, flag)
%% Animation and Solution of Double Pendulum Motion

% * the length of the first rod $L_1$
% * the length of the second rod $L_2$
% * the mass of the first bob $m_1$
% * the mass of the second bob $m_2$
% * the gravitational constant $g$

syms theta_1(t) theta_2(t) L_1 L_2 m_1 m_2 g
%% 
% Define the displacements of the double pendulum in Cartesian coordinates.

x_1 = L_1*sin(theta_1);
y_1 = -L_1*cos(theta_1);
x_2 = x_1 + L_2*sin(theta_2);
y_2 = y_1 - L_2*cos(theta_2);

% Find the velocities by differentiating the displacements with respect to time 

vx_1 = diff(x_1);
vy_1 = diff(y_1);
vx_2 = diff(x_2);
vy_2 = diff(y_2);
 
% Find the accelerations by differentiating the velocities with respect to time.

ax_1 = diff(vx_1);
ay_1 = diff(vy_1);
ax_2 = diff(vx_2);
ay_2 = diff(vy_2);
%% Step 2: Define Equations of Motion
% Define the equations of motion based on Newton's laws.
% 
% First, specify the tension of the first rod as $T_1$, and the tension of the 
% second rod $T_2$.

syms T_1 T_2
%% 
% Next, construct free-body diagrams of the forces that act on both masses.
% 
% 
% 
% Evaluate the forces acting on $m_1$. Define the equations of motion of the 
% first bob by balancing the horizontal and vertical force components. Specify 
% these two equations as symbolic equations |eqx_1| and |eqy_1|.

eqx_1 = m_1*ax_1(t) == -T_1*sin(theta_1(t)) + T_2*sin(theta_2(t));
eqy_1 = m_1*ay_1(t) == T_1*cos(theta_1(t)) - T_2*cos(theta_2(t)) - m_1*g;
% 
% Evaluate the forces acting on $m_2$. Define the equations of motion of the 
% second bob by balancing the horizontal and vertical force components. Specify 
% these two equations as symbolic equations |eqx_2| and |eqy_2|.

eqx_2 = m_2*ax_2(t) == -T_2*sin(theta_2(t));
eqy_2 = m_2*ay_2(t) == T_2*cos(theta_2(t)) - m_2*g;
%% Step 3: Evaluate Forces and Reduce System Equations
% Four equations of motion describe the kinematics of the double pendulum. Evaluate 
% the forces acting on the rods and reduce the set of four equations to two equations.
% 
% The equations of motion have four unknowns: $\theta_1$, $\theta_2$, $T_1$, 
% and $T_2$. Evaluate the two unknowns $T_1$ and $T_2$ from |eqx_1| and |eqy_1|. 
% Use <docid:symbolic_ug#buezrr6 solve> function to find $T_1$ and $T_2$.

Tension = solve([eqx_1 eqy_1],[T_1 T_2]);
%% 
% Substitute the solutions for $T_1$ and $T_2$ into |eqx_2| and |eqy_2|.

eqRed_1 = subs(eqx_2,[T_1 T_2],[Tension.T_1 Tension.T_2]);
eqRed_2 = subs(eqy_2,[T_1 T_2],[Tension.T_1 Tension.T_2]);
%% 
% The two reduced equations fully describe the pendulum motion.
%% Step 4: Solve System Equations
% Solve the system equations to describe the pendulum motion.
% 
% First, define the values for the masses in $\textrm{kg}$, the rod lengths 
% in $\mathrm{m}$, and the gravity in $\mathrm{m}/{\mathrm{s}}^2$ (SI units). 
% Substitute these values into the two reduced equations.

L_1 = 1;
L_2 = 1.5;
m_1 = 2;
m_2 = 1;
g = 9.8;
eqn_1 = subs(eqRed_1);
eqn_2 = subs(eqRed_2);
%% 
% The two equations are nonlinear second-order differential equations. To solve 
% these equations, convert them to first-order differential equations.

[V,S] = odeToVectorField(eqn_1,eqn_2);
%% 
% The elements of the vector |V| represent the first-order differential equations 
% that are equal to the time derivative of the elements of |S|. The elements of 
% |S| are the state variables $\theta_2$, $d\theta_2 /\textrm{dt}$, $\theta_1$, 
% and $d\theta_1 /\textrm{dt}$. The state variables describe the angular displacements 
% and velocities of the double pendulum.

%% 
% Next, convert the first order-differential equations to a MATLAB function 
% with the handle |M|.

M = matlabFunction(V,'vars',{'t','Y'});
%% 
% Define the initial conditions of the state variables as |[pi/4 0 pi/6 0]|. 
% Use the <docid:matlab_ref#bu00_4l ode45> function to solve for the state variables. 
% The solutions are a function of time within the interval |[0 10]|.

[T,Y] = ode45(M,tspan,initCond);
if flag
    %% 
    % Plot the solutions of the state variables.
    figure;
    plot(T,Y)
    legend('\theta_2','d\theta_2/dt','\theta_1','d\theta_1/dt')
    title('Solutions of State Variables')
    xlabel('Time (s)')
    ylabel('Solutions (rad or rad/s)')
end
end
% %% Step 5: Create Animation of Oscillating Double Pendulum
% % Create the animation of the oscillating double pendulum.
% % 
% % First, create four functions that use deval to 
% % evaluate the coordinates of both pendulums from the solutions |sols|.
% 
% x_1 = @(t) L_1*sin(deval(sols,t,3));
% y_1 = @(t) -L_1*cos(deval(sols,t,3));
% x_2 = @(t) L_1*sin(deval(sols,t,3))+L_2*sin(deval(sols,t,1));
% y_2 = @(t) -L_1*cos(deval(sols,t,3))-L_2*cos(deval(sols,t,1));
% %% 
% % Next, create a stop-motion animation object of the first pendulum bob by using 
% % the fanimator function. 
% % By default, |fanimator| creates an animation object with 10 generated frames 
% % per unit time within the range of |t| from 0 to 10. Plot the coordinates by 
% % using the |plot| function. Set the _x_-axis and _y_-axis to be equal length.
% figure;
% fanimator(@(t) plot(x_1(t),y_1(t),'ro','MarkerSize',m_1*10,'MarkerFaceColor','r'),...
%         'AnimationRange', [0, tSim]);
% axis equal;
% %% 
% % Next, add the animation objects of the first rigid rod, the second pendulum 
% % bob, and the second rigid rod.
% 
% hold on;
% fanimator(@(t) plot([0 x_1(t)],[0 y_1(t)],'r-'), 'AnimationRange', [0, tSim]);
% fanimator(@(t) plot(x_2(t),y_2(t),'go','MarkerSize',m_2*10,'MarkerFaceColor','g'),...
%         'AnimationRange', [0, tSim]);
% fanimator(@(t) plot([x_1(t) x_2(t)],[y_1(t) y_2(t)],'g-'),...
%         'AnimationRange', [0, tSim]);
% fanimator(@(t) text(-0.3,0.3,"Timer: "+num2str(t,4)), 'AnimationRange', [0, tSim]);
% hold off;
% playAnimation;