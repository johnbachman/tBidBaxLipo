
nbd_mutants = ['c3', 'c62', 'c120', 'c122', 'c126']

def linear_pathway(num_steps=5, num_phases=1):
    ode_output = "function out = nbd_odes(t, input, param)\n\n"
    ode_species = ''
    odes = ''
    ode_params = ''
    ode_param_counter = 1


    nbd_output = 'function [t nbd_out] = nbd_run(nbd_params, time)\n\n'
    nbd_initvals = 'initvals = ['
    #nbd_paramvals = ''
    nbd_params = ''
    nbd_param_counter = 1
    nbd_eqns = ''


    # Parameters: There will be two parameters (forward and reverse) for each
    # step (so for two steps, four parameters, etc.)
    # In addition, for each mutant, there will be a parameter for each of the
    # "phases" it is presumed to depend on. 

    # 1. Create the mapping between names and parameters for the ODE file
    for step_num in range(1, num_steps+1):
        ode_params += 'kf%d = param(%d);\n' % (step_num, ode_param_counter)
        ode_params += 'kr%d = param(%d);\n' % (step_num, ode_param_counter + 1)
        nbd_params += 'param(%d) = nbd_params(%d);\n' % (ode_param_counter, nbd_param_counter)
        nbd_params += 'param(%d) = nbd_params(%d);\n' % (ode_param_counter + 1, nbd_param_counter + 1)
        ode_param_counter += 2
        nbd_param_counter += 2
    # 2. Create the mapping between names and parameters for the NBD file
    for i, nbd_mutant in enumerate(nbd_mutants):
        for j in range(1, num_phases+1):
            nbd_params += 'k_%s_%d = nbd_params(%d);\n' % (nbd_mutant, j, nbd_param_counter)
            nbd_param_counter += 1

    #for i in range(1, param_counter):
    #    paramvals += 'param(%d) = 1;\n' % i

    # ODEs for state S[i]
    num_states = num_steps + 1
    for state_num in range(1, num_states+1):
        ode_species += 'S%d = input(%d);\n' % (state_num, state_num)
        # These equations will have terms arising from the other states S[i]
        # in a purely linear order
        odes += '%% S%d\n' % state_num        
        # Equations differ for first, last, and intermediate states in ther series
        if (state_num == 1):
            odes += 'out(%d, 1) = -(kf%d * S%d) + (kr%d * S%d);\n' % \
                    (state_num, state_num, state_num, state_num, state_num + 1)
            nbd_initvals += '1 '
        elif (state_num == num_states):
            odes += 'out(%d, 1) = (kf%d * S%d) - (kr%d * S%d);\n' % \
                    (state_num, state_num-1, state_num-1, state_num-1, state_num)
            nbd_initvals += '0 '
        else:
            odes += 'out(%d, 1) = -((kf%d + kr%d) * S%d) + (kf%d * S%d) + (kr%d * S%d);\n' % \
                    (state_num,
                     state_num, state_num-1, state_num,
                     state_num-1, state_num-1,
                     state_num, state_num + 1)
            nbd_initvals += '0 '

    # Eqns for NBD signals
    #== K1 ==
    #3c:   1.381e-2,  SD 9.050e-4
    #122c: 1.159e-2,  SD 1.140e-3
    #62c:  6.356e-3,  SD 1.567e-3 (nearly worthless)
    #126c: 5.875e-3
    #120c: 1.981e-3,  SD 1.233e-4
    # Note that these are 1-indexed, not 0-indexed!

    signal_states = {'c3':[2], 'c122':[3], 'c62':[4], 'c126':[5], 'c120':[6]}

    for i, nbd_mutant in enumerate(nbd_mutants):
        #species_index = num_states + 1 + i
        #species += '%s = input(%d);\n' % (nbd_mutant, species_index)
        #odes += '%% %s\n' % nbd_mutant
        # These equations will have terms arising only from the states S[i], not
        # from any of the other NBD equations

        # Choose which states will contribute to the signal of this mutant
        #signal_states = [2, 3] # FIXME This should be generated systematically
        # FIXME This should be generated systematically
        
        start_vals = [1.0, 1.0656, 1.2964, 0.9884, 1.0] # FIXME this should not be hard-coded!

        nbd_eqns += '%% %s\n' % (nbd_mutant)
        nbd_eqns += 'nbd_out(:,%d) = ' % (i+1)
        last_signal_state = str(start_vals[i])
        for j in range(1, num_phases+1):
            # the state for this mutant
            signal_states[nbd_mutant].sort()
            jth_signal_state = signal_states[nbd_mutant][j-1]
            # If this is the first phase for this mutant, i.e., j = 1, then all terms less
            # than this one should have the preceding term;
            if (j == 1):
                for k in range(1,jth_signal_state):
                    nbd_eqns += '(%s * x(:,%d)) + ' % (last_signal_state, k)

            # The signal value for this signal state
            signal_state_param = 'k_%s_%d' % (nbd_mutant, j)
            nbd_eqns += '(%s * x(:,%d))' % (signal_state_param, jth_signal_state)
            if (jth_signal_state == num_steps+1):
                nbd_eqns += ';\n'
            else:
                nbd_eqns += ' + '

            # and if this is the last signal phase for this mutant, then all
            # subsequent ones should have this weight ascribed to them
            if (j == num_phases):
                for k in range(jth_signal_state+1, num_steps+2):
                    nbd_eqns += '(%s * x(:,%d))' % (signal_state_param, k)
                    if (k == num_steps+1):
                        nbd_eqns += ';\n'
                    else:
                        nbd_eqns += ' + '

            # Save the last signal state
            last_signal_state = signal_state_param

    nbd_initvals += '];\n'

    ode_output += ode_params + '\n' + ode_species + '\n' + odes + '\nend\n'

    nbd_output += nbd_params + '\n'
    nbd_output += nbd_initvals + '\n'
    nbd_output += '% Set the time length of the experiment\n'
    nbd_output += '% tspan = [0:2:3590];\n\n'
    nbd_output += "% Accuracy of our integrator. Don't worry about this line.\n"
    nbd_output += "options=odeset('AbsTol', 1e-15, 'RelTol', 2.22045e-14);\n\n"
    nbd_output += "[t x] = ode15s(@nbd_odes, time, initvals, options, param);\n\n"
    nbd_output += nbd_eqns
    nbd_output += '\n\nend'

    print "NBD OUTPUT ========================\n"
    print nbd_output

    print "ODE OUTPUT ========================\n"
    print ode_output
    print '\n\n'
    #print paramvals
    ode_output_file = open('nbd_odes.m', 'w')
    ode_output_file.write(ode_output)
    ode_output_file.close()

    nbd_output_file = open('nbd_run.m', 'w')
    nbd_output_file.write(nbd_output)
    nbd_output_file.close()


linear_pathway()

# ODE for each of the residues

# ODE for each of the states

