def cv_51_minute_plots():
    '''
    latest version
    '''


    import mesa_reader as mr

    paths = ['history_2p63days.data', 'history_2p69days.data', 'history_2p75days.data']
    colors = ['r', 'k', 'b']
    p_inits = [2.63, 2.69, 2.75]
    t_ms = 7583974667 # age at which a 1.1 Msun star passes the Schoenberg-Chandrasekhar limit and leaves MS
    
    nrow, ncol = 5, 2
    f, axx = plt.subplots(nrow, ncol, figsize = (12, 17))
    plt.subplots_adjust(hspace = 0, wspace=0.2)
    ax = np.concatenate(axx)
    for i in range(nrow):
        for j in range(ncol):
            axx[i, j].set_xlim(0, 24)
            if i != nrow -1:
                axx[i, j].set_xticklabels([])

    t_over_tms_at_rlof = []
    for i, path in enumerate(paths):
        h1 = mr.MesaData(path)
        rlof = (10**h1.log_R/h1.rl_1 > 0.99)
        t_over_tms_at_rlof.append(np.min(h1.star_age[rlof])/t_ms)
    
        label_, = ax[0].plot(24*h1.data('period_days'), h1.data('star_mass'), color= colors[i],label = r'$P_{\rm init}=%.2f\,\rm days$' % p_inits[i])
        ax[1].plot(24*h1.data('period_days'), 10**h1.data('log_R'), color= colors[i])
        ax[1].plot(24*h1.data('period_days'), h1.data('rl_1'), ls=':', color= colors[i])
        ax[2].plot(24*h1.data('period_days'), 10**h1.data('log_Teff'), color= colors[i], label = r'$t_{\rm RLOF}/t_{\rm MS}= %.2f$' % t_over_tms_at_rlof[i])
        ax[3].plot(24*h1.data('period_days'), h1.data('total_mass_h1'), color= colors[i])
        ax[4].plot(24*h1.data('period_days'), h1.data('star_age')/1e9, color= colors[i])
        
        mdot_smooth = running_median_insort(h1.data('lg_mstar_dot_1'), 10)
    
        t_accretor = 1.7e4*(10**mdot_smooth/1e-10)**(1/4)*(0.55/0.9) # assuming Teff set by compressional heating
    
        ax[5].plot(24*h1.data('period_days'), mdot_smooth, color= colors[i])
        ax[6].plot(24*h1.data('period_days'), t_accretor, color= colors[i])
     
        # MB is suppressed when q_conv < 0.02
        jdot_mb_pred = -6.8e34*h1.star_mass * (10**h1.log_R)**3 * (h1.period_days)**-3
        qconv = h1.conv_mx1_top - h1.conv_mx1_bot 
        supress = np.exp(1-0.02/qconv)
        supress[qconv > 0.02] = 1
    
        ax[7].plot(24*h1.data('period_days'), np.log10(-h1.data('jdot_mb')), colors[i])
        ax[7].plot(24*h1.data('period_days'), np.log10(-h1.data('jdot_gr')), colors[i], ls = ':')        
        ax[8].plot(24*h1.data('period_days'), h1.surface_he4, colors[i])
        solar_n_to_c = h1.surface_n14[0]/h1.surface_c12[0]
        ax[9].plot(24*h1.data('period_days'), np.log10((h1.surface_n14/h1.surface_c12)/solar_n_to_c), colors[i])        
    
    ax[1].plot([], [], ':', color = '0.5', label = r'$R_{\rm Roche\,\,lobe}$')
    ax[1].set_ylim(0, 2.3)
    ax[0].set_ylim(0, 1.33)
    ax[3].set_ylim(-0.01, 0.63)

    ax[2].set_ylim(3.2e3, 7.3e3)
    ax[5].set_ylim(-12, -6.5)
    ax[6].set_ylim(1000, 58000)
    ax[7].set_ylim(30, 37.5)
    ax[4].set_ylim(6.5, 8.5)

    ax[7].plot([], [], ':', color = '0.5', label = r'$\dot{J}_{\rm GR}$')
    ax[7].plot([], [], '-', color = '0.5', label =  r'$\dot{J}_{\rm MB}$')
    ax[7].legend(loc = 'lower right', frameon= False, fontsize=16)
 
    ax[1].legend(loc = 'upper left', frameon = False, fontsize=16)
    ax[2].legend(loc = 'upper right', frameon = False,  fontsize=10)

    axx[nrow-1, 0].set_xlabel(r'$P_{\rm orb}\,[\rm hr]$', fontsize=20)
    axx[nrow-1, 1].set_xlabel(r'$P_{\rm orb}\,[\rm hr]$', fontsize=20)

    ax[0].set_ylabel(r'$M_{\rm donor}\,[M_{\odot}]$', fontsize=20)
    ax[1].set_ylabel(r'$R_{\rm donor}\,[R_{\odot}]$', fontsize=20)
    ax[2].set_ylabel(r'$T_{\rm eff,\,donor}\,[\rm K]$', fontsize=20)
    ax[3].set_ylabel(r'$M_{\rm H}\,\,[M_{\odot}]$', fontsize=20)
    ax[4].set_ylabel(r'$\rm time \,\,[\rm Gyr]$', fontsize=20)
    ax[5].set_ylabel(r'$\log(\dot{M}/(M_{\odot}\,\rm yr^{-1}))$', fontsize=20)
    ax[6].set_ylabel(r'$T_{\rm eff,\,accretor}\,[\rm K]$', fontsize=20)
    ax[7].set_ylabel(r'$\log\left(-\dot{J}/{\rm \left(dyn\,cm\right)}\right)$', fontsize=20)
    ax[8].set_ylabel(r'$Y_{\rm He,\,donor\,surface}$', fontsize=20)
    ax[9].set_ylabel(r'$\rm [N/C]_{\rm donor\,surface}$', fontsize=20)

    for i in range(nrow):
        for j in range(ncol):
            yticks = axx[i, j].get_ylim()
            axx[i, j].plot(2*[51/60.], yticks, 'k--', lw=0.5)
            axx[i, j].set_ylim(yticks)
    ax[0].plot([51/60], [0.1184], 'c*', ms=10)
    ax[1].plot([51/60], [0.1017], 'c*', ms=10)
    ax[2].plot([51/60], [6000], 'c*', ms=10)
    ax[6].plot([51/60], [12600], 'c*', ms=10)

def running_median(seq, window_size = 30):
    '''
    this is just a fast way of calculating a running median
    '''
    from collections import deque
    from bisect import insort, bisect_left
    from itertools import islice
    
    seq = iter(seq)
    d = deque()
    s = []
    result = []
    for item in islice(seq, window_size):
        d.append(item)
        insort(s, item)
        result.append(s[len(d)//2])
        m = window_size // 2
    for item in seq:
        old = d.popleft()
        d.append(item)
        del s[bisect_left(s, old)]
        insort(s, item)
        result.append(s[m])
    return np.array(result)