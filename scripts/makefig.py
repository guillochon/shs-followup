import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / 'paper'

G=6.674e-8; c=2.998e10; Msun=1.989e33; pc=3.086e18; km=1e5
def rg(M1): return G*M1*Msun/c**2
def v_inf(M1,q12,rp):
    M2=M1/q12; rH=rp*(M2/(3*M1))**(1/3.)
    return 1.3*q12**(1/6.)*np.sqrt(G*M2*Msun/rH)

plt.rcParams.update({'font.size':9,'font.family':'serif','mathtext.fontset':'dejavuserif',
                     'axes.linewidth':0.8,'xtick.direction':'in','ytick.direction':'in',
                     'xtick.top':True,'ytick.right':True,'legend.frameon':False,
                     'pdf.fonttype':42,'ps.fonttype':42})

#=========================== FIGURE 1 ===================================
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(7.1,3.0))

# ---------------- Panel 1 ----------------
M1,q12=1e9,10.
e=1.-np.logspace(0,-4.5,4000)   # dense near e=1 so every curve reaches c/2
cols=['#1b4965','#3d8ea8','#7fb3c8']
vlim=0.5*c/km
for a_pc,col,ls in zip([0.05,5e-3,1e-3],cols,['-','-','-']):
    rp=a_pc*pc*(1-e)
    v=v_inf(M1,q12,rp)/km
    v=np.where(rp/rg(M1)>4.,v,np.nan)          # IBCO floor: r_p > 4 r_g
    v=np.minimum(v,vlim)                        # GL15 speed limit c/2
    ax1.plot(e,v,ls,color=col,lw=1.7,label=r'$a=%g$ pc'%a_pc)

# accessible band, and the swallowed region above c/2 (shaded grey)
ax1.axhspan(2.4e4,vlim,color='#e8b04b',alpha=0.16,zorder=0)
ax1.axhspan(vlim,3.0e5,color='0.55',alpha=0.20,zorder=0)
ax1.axhline(vlim,color='0.45',lw=0.9)
ax1.text(0.20,1.13e5,'accessible only to\nmerging MBH binaries',
         fontsize=6.0,color='#8a6410',ha='center',va='center',linespacing=1.5)
ax1.text(0.03,2.15e5,'stars swallowed\n($r_{p}<r_{\\rm IBCO}$, $v_\\infty>c/2$)',
         fontsize=6.2,color='0.25',ha='left',va='center',linespacing=1.5)

# Hills-mechanism ceilings; labels sit clear of both the lines and the curves
ax1.axhline(1.1e4,color='0.45',ls=':',lw=0.9)
ax1.axhline(2.4e4,color='0.45',ls='--',lw=0.9)
ax1.text(0.03,7.5e3,'star$+$MBH limit',fontsize=6.2,color='0.35',ha='left',va='center')
ax1.text(0.03,1.75e4,'SBH$+$star limit',fontsize=6.2,color='0.35',ha='left',va='center')

# Hardening-phase e from the cited merger channels: prograde CBD attractor
# (Zrake/Siwek/Muñoz, 0.45–0.5) plus stellar scattering and retrograde disks
# (Iwasawa/Khan/Tiede), which extend toward e → 1.
ax1.axvspan(0.45,1.00,color='#c1443c',alpha=0.18,zorder=0)
ax1.axvspan(0.45,0.50,color='#c1443c',alpha=0.32,zorder=0)
ax1.text(0.475,1.35e3,'CBD\n$e_{\\rm eq}$',fontsize=6.2,color='#8f2f28',ha='center',va='bottom')
ax1.text(0.78,1.35e3,'stellar /\nretrograde',fontsize=6.2,color='#8f2f28',ha='center',va='bottom')
ax1.set_yscale('log'); ax1.set_xlim(0,1.0); ax1.set_ylim(1e3,3.0e5)
ax1.set_xlabel(r'binary eccentricity $e_{12}$')
ax1.set_ylabel(r'$v_\infty$ (km s$^{-1}$)')
ax1.legend(loc='lower left',fontsize=7,handlelength=1.5,
           bbox_to_anchor=(0.02,0.02),borderaxespad=0)
ax1.set_title('(a) ejection speed vs. eccentricity',fontsize=8.5,pad=4)

# ---------------- Panel 2 ----------------
wb=dict(facecolor='white',edgecolor='none',pad=0.8)
mlim=np.linspace(18.8,29,200)
for M,lab,col in [(-6.0,'bright TP-AGB','#1b4965'),(-4.5,'TRGB','#c1443c')]:
    d=10**((mlim-M)/5.+1)/1e6
    ax2.plot(mlim,d,color=col,lw=1.7,label=lab)
ax2.axhline(16,color='0.45',ls='--',lw=0.9)
ax2.text(19.15,19.5,r'$r_{\rm Virgo}=16$ Mpc',fontsize=6.4,color='0.35',va='bottom',bbox=wb)
for m,lab,ty in [(19.5,'Euclid\ngrism',0.105),(21,'HLWAS\ngrism',0.022),
                 (24,'Euclid\nphot.',0.022),(26.2,'HLWAS\nwide',0.022),
                 (26.5,'HLWAS\nmedium',0.105),(27.5,'HLWAS\ndeep',0.022)]:
    ax2.axvline(m,color='0.7',lw=0.8,ls=':')
    ax2.text(m,ty,lab,fontsize=6.0,ha='center',va='bottom',color='0.4',bbox=wb)
ax2.set_yscale('log'); ax2.set_xlim(18.8,29); ax2.set_ylim(1.5e-2,1e2)
ax2.set_xlabel(r'limiting magnitude $m_{\rm AB}$ (near-IR)')
ax2.set_ylabel(r'maximum distance (Mpc)')
ax2.legend(loc='upper left',fontsize=7,handlelength=1.5)
ax2.set_title(r'(b) distance reach',fontsize=8.5,pad=4)

fig.tight_layout(pad=0.5,w_pad=1.6)
fig.savefig(PAPER/'f1.pdf',bbox_inches='tight')
print("wrote", PAPER/'f1.pdf')

#=========================== FIGURE 2 ===================================
# Expected detections.  Fiducial is the unboosted GL15 density (one SHS per
# 100 deg^2 at the H<24 Euclid/WFIRST-class limit; those papers quote AB).
# That density already assumes e≳0.5 plunging experiments, so the 1.7 from
# P(v) is not applied on top.  Euclidean growth is the primary prediction.
# The shaded band is the eccentricity systematic only (e=0 to e=0.9).
ECC_LO=1.7                    # e=0 sits 1.7 below the (already e≳0.5) fiducial
ECC_HI=6.                     # e=0.9 count boost relative to circular; retrograde
S24=0.01                      # deg^-2 at m_AB = 24, unboosted GL15
A_REF=2356.                   # deg^2, HLWAS medium / grism (STScI 2026 Aug)
RV=16.                        # Mpc
M_TRGB, M_TPAGB = -4.5, -6.0
F_TRGB = F_TPAGB = 0.5        # GL15 detectable split at H=24
# Ivezic+08 Table 2 is in r; cool-giant colours onto F158 (AB).
R_MINUS_F158_TRGB = 1.5
R_MINUS_F158_TPAGB = 2.5
# ETC-scaled HLWAS medium-tier stack (4 passes x 2 dithers x 190 s):
# 21.3 (5σ/pix, 1 hr, 2x zodi) -> 20.8; ~30% lower HLWAS zodi -> 20.9.
# Marker at the ≈21 adopted in the text; bar spans the 1-hr published range.
GRISM_M, GRISM_M_LO, GRISM_M_HI = 21.0, 20.5, 21.3
VINF = 1e4                    # km s^-1, SHS threshold
A_WIDE, A_MED, A_DEEP = 2525., 2356., 19.2
A_HLWAS = A_WIDE + A_MED      # 4881 deg^2
A_EUCLID = 14000.
M_DEEP = 27.5                 # F158 deep-tier 5σ (STScI HLWAS table)

def dist(m, M):
    return 10**((np.asarray(m, dtype=float)-M)/5. + 1)/1e6

def mag(d_mpc, M):
    return M + 5*np.log10(np.asarray(d_mpc, dtype=float)*1e6) - 5

def sigma_mu_ivezic(r):
    """Ivezic+08 Table 2: 10-yr LSST proper-motion 1σ (mas yr^-1)."""
    r = np.asarray(r, dtype=float)
    r_tab = np.array([21., 22., 23., 24.])
    s_tab = np.array([0.2, 0.3, 0.5, 1.0])
    out = np.empty_like(r)
    bright = r <= 21
    out[bright] = 0.2
    mid = (r > 21) & (r <= 24)
    out[mid] = np.interp(r[mid], r_tab, s_tab)
    faint = r > 24
    out[faint] = 1.0 * 10**(0.3*(r[faint]-24))
    return float(out) if out.ndim == 0 else out

def rubin_dmax(M, color, nsig=5., v=VINF):
    """Farthest d (Mpc) at which mu(v,d) >= nsig * sigma_mu(r), r = F158+color."""
    ds = np.logspace(-2, 1.5, 4000)
    mF = mag(ds, M)
    r = mF + color
    mu = v / (4.74 * ds * 1000.)
    ok = (mu >= nsig * sigma_mu_ivezic(r)) & (mF < 26)
    return float(ds[ok][-1]) if np.any(ok) else 0.

def sigma_one(m, M, f, rmax=None):
    """one class: n from f*S24 at that tracer's m=24 horizon; N ∝ d^3"""
    d24 = dist(24., M)
    d = dist(m, M)
    if rmax is not None:
        d = np.minimum(d, rmax)
    return f*S24*(d/d24)**3

def sigma_euc(m):
    return S24*10**(0.6*(np.asarray(m, dtype=float)-24.))

def sigma_virgo(m):
    return sigma_one(m, M_TRGB, F_TRGB, RV) + sigma_one(m, M_TPAGB, F_TPAGB, RV)

def n_distance_limited(A, dT, dP):
    d24_T, d24_P = dist(24., M_TRGB), dist(24., M_TPAGB)
    return A*(F_TRGB*S24*(dT/d24_T)**3 + F_TPAGB*S24*(dP/d24_P)**3)

def fmt_n(n):
    n=float(np.squeeze(n))
    if n>=100: return '%d'%(int(n/10+0.5)*10)
    if n>=10:  return '%d'%int(round(n))
    if n>=1:   return '%d'%int(round(n))
    return '%.1f'%n

d_pm_trgb = rubin_dmax(M_TRGB, R_MINUS_F158_TRGB)
d_pm_tp   = rubin_dmax(M_TPAGB, R_MINUS_F158_TPAGB)
n_rubin_hlwas = n_distance_limited(A_HLWAS, d_pm_trgb, d_pm_tp)
n_rubin_lsst  = n_distance_limited(18000., d_pm_trgb, d_pm_tp)
m_pm = mag(d_pm_trgb, M_TRGB)   # characteristic F158 of the 5σ TRGB cut in r
# HLTDS: 18 deg^2, 5-yr Roman-only baseline, σ_μ ≈ 0.1 mas/yr in F158
SIG_HLTDS = 0.1
d_hltds = VINF / (4.74 * 1000. * 5. * SIG_HLTDS)
n_hltds = n_distance_limited(18., d_hltds, d_hltds)
print('Rubin 5-sigma dmax TRGB/TP-AGB = %.2f / %.2f Mpc'%(d_pm_trgb, d_pm_tp))
print('Rubin N HLWAS / LSST = %.2f / %.2f'%(n_rubin_hlwas, n_rubin_lsst))
print('Rubin characteristic m_F158 = %.2f'%m_pm)
print('HLTDS dmax / N = %.2f Mpc / %.3f'%(d_hltds, n_hltds))
print('Euclid phot / wide / med / deep = %.1f / %.1f / %.1f / %.1f'%(
    A_EUCLID*sigma_euc(24.0), A_WIDE*sigma_euc(26.2),
    A_MED*sigma_euc(26.5), A_DEEP*sigma_euc(M_DEEP)))
print('Grism continuum Euclid / Roman = %.3f / %.3f'%(
    A_EUCLID*sigma_euc(19.5), A_MED*sigma_euc(GRISM_M)))

# Envelope about the unboosted GL15 fiducial: e=0 is /1.7, e=0.9 is x6
FAC_LO, FAC_HI = ECC_LO, ECC_HI

fig2,ax=plt.subplots(figsize=(3.5,3.35))
m=np.linspace(19,29,400)
y_euc=A_REF*sigma_euc(m)
y_cap=A_REF*sigma_virgo(m)

ax.fill_between(m, y_euc/FAC_LO, y_euc*FAC_HI, color='#c1443c', alpha=0.16, lw=0,
                zorder=0, label=r'eccentricity only ($e=0$--$0.9$)')
ax.plot(m, y_euc, color='#1b2838', lw=1.9, label='Euclidean (primary)')
ax.plot(m, y_cap, color='0.45', ls='--', lw=1.4,
        label=r'Virgo floor ($r\leq16$ Mpc)')

# survey markers sit at each survey's own (m, N); N is Euclidean except Rubin (PM)
surveys=[
    # mm, N, name, area, marker, kind
    # Published 3.5σ continuum H=19.5 ± 0.2 (Gabarra+23); two-sided like Roman.
    (19.5, float(A_EUCLID*sigma_euc(19.5)), 'Euclid grism', '14,000 deg$^2$', 'D'),
    (GRISM_M, float(A_MED*sigma_euc(GRISM_M)), 'Roman grism', '2,350 deg$^2$', 'o'),
    (m_pm, n_rubin_hlwas, 'Rubin PM', 'HLWAS overlap', 's'),
    (24.0, float(A_EUCLID*sigma_euc(24.0)), 'Euclid phot.', '14,000 deg$^2$', '^'),
    (26.2, float(A_WIDE*sigma_euc(26.2)), 'HLWAS wide', '2,500 deg$^2$', 'v'),
    (26.5, float(A_MED*sigma_euc(26.5)), 'HLWAS medium', '2,350 deg$^2$', 'v'),
    (M_DEEP, float(A_DEEP*sigma_euc(M_DEEP)), 'HLWAS deep', '19.2 deg$^2$', 'v'),
]
for mm,y,name,area,mk in surveys:
    ylo,yhi=y/FAC_LO, y*FAC_HI
    xerr = None
    if name=='Roman grism':
        mfc,ms='#22c55e',6.0
        xerr = np.array([[mm-GRISM_M_LO],[GRISM_M_HI-mm]])
        xy,ha,va=(21.45, ylo*0.72),'left','top'
    elif name=='Euclid grism':
        mfc,ms='#86efac',5.4
        xerr = np.array([[0.2],[0.2]])  # published 19.5 ± 0.2
        xy,ha,va=(19.08, yhi*1.15),'left','bottom'
    elif mk=='s':
        mfc,ms='#0e7c7b',5.0
        xy,ha,va=(mm+0.18, yhi*1.15),'left','bottom'
    elif name=='Euclid phot.':
        mfc,ms='#a7f3d0',5.2
        xy,ha,va=(mm-0.18, y),'right','center'
    elif name=='HLWAS wide':
        mfc,ms='w',5.2
        xy,ha,va=(mm-0.18, y),'right','center'
    elif name=='HLWAS medium':
        mfc,ms='w',5.2
        xy,ha,va=(mm-0.10, yhi*1.12),'right','bottom'
    else:
        mfc,ms='w',5.2
        xy,ha,va=(mm-0.18, y*0.55),'right','top'
    ax.errorbar(mm,y,xerr=xerr,yerr=np.array([[y-ylo],[yhi-y]]),
                fmt=mk,ms=ms,mfc=mfc,mec='0.2',mew=1.1,
                ecolor='0.3',elinewidth=0.8,capsize=2.0,capthick=0.8,zorder=5)
    lab='%s\n%s\n%s SHS'%(name,area,fmt_n(y))
    ax.text(xy[0],xy[1],lab,ha=ha,va=va,zorder=6,
            fontsize=5.2,color='0.15',linespacing=0.88,clip_on=True,
            bbox=dict(facecolor='white',edgecolor='none',pad=0.12,alpha=0.88))
    print('%s: m=%.2f  N=%.2f  (%.2f--%.2f ecc)'%(name, mm, y, ylo, yhi))

ax.set_yscale('log'); ax.set_xlim(19,29); ax.set_ylim(3e-3,1e5)
ax.set_xlabel(r'limiting magnitude $m_{\rm AB}$ (near-IR)')
ax.set_ylabel('expected number of SHS')
ax.set_title('expected survey yield',fontsize=8.5,pad=4)
leg=ax.legend(loc='lower right',fontsize=5.4,handlelength=2.4,labelspacing=0.20,
              borderpad=0.15,handletextpad=0.35,frameon=False)
leg.get_frame().set_facecolor('none')
leg.get_frame().set_edgecolor('none')
fig2.subplots_adjust(left=0.18,right=0.98,top=0.90,bottom=0.14)

fig2.savefig(PAPER/'f2.pdf',bbox_inches='tight')
print("wrote", PAPER/'f2.pdf')
