import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import date
from skyfield.api import load
from astropy.time import Time

vernal_zps = [59658.64791667, 60023.89166667, 60389.12916667, 60754.37569444, 61119.61527778]
autumnal_zps = [59845.04444444, 60210.28472222, 60575.53055556, 60940.76319444, 61306.00347222]

planets = load('./de440s.bsp')
earth = planets['earth']
moon = planets['moon']
ts = load.timescale()

linestyles = {0:'-', 1:'-.', 2:':', 3:'--'}

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman'
plt.rcParams['mathtext.bf'] = 'Times New Roman'
plt.rcParams.update({'font.size': 16})


today = str(date.today()) + 'T00:00.0'
t = Time(today, format = 'isot', scale = 'utc')
today_mjd = t.mjd

def sin(angle):
    return np.sin(np.deg2rad(angle))
def cos(angle):
    return np.cos(np.deg2rad(angle))

def arcsin(value):
    return np.rad2deg(np.arcsin(value))
def arccos(value):
    return np.rad2deg(np.arccos(value))

class observatory():
    def __init__(self):
        self.latitude, self.longitude = 0.0, 0.0
        self.id = ''
        self.name = ''
        self.country = ''
        self.elevation = 0.0
        self.light_pollution = ''
        self.limiting_magnitude = 0.0
        
def convert_latitude(angle):
    split = angle.split()
    if split[-1] == 'S':
        factor = -1
    else:
        factor = 1
        
    return factor * (float(split[0][:-1]) + float(split[1][:-1])/60 + float(split[2][:-1])/3600)

def convert_longitude(angle):
    split = angle.split()
    if split[-1] == 'W':
        factor = -1
        offset = 360
    else:
        factor = 1
        offset = 0
        
    return offset + factor * (float(split[0][:-1]) + float(split[1][:-1])/60 + float(split[2][:-1])/3600)

def retrieve_observatory(iden):
    data = pd.read_csv('./observatories.csv', delimiter = '\t')
    if iden not in list(data['ID']):
        print('The observatory ID \'' + iden + '\' was not found in the database. Known IDs are ' + str(list(data['ID'])) + '.')
        print('\nAlternatively you may define a new observatory.')
        return -1
    else:
        obs = observatory()
        for i in range(len(data['ID'])):
            if iden == data['ID'][i]:
                obs.id = iden
                obs.name = data['Name'][i]
                obs.latitude = convert_latitude(data['Latitude'][i])
                obs.longitude = convert_longitude(data['Longitude'][i])
                obs.country = data['Country'][i]
                obs.light_pollution = data['Light pollution'][i]
                obs.limiting_magnitude = data['Limiting magnitude'][i]
                obs.elevation = data['Elevation'][i]
                continue
        return obs

def convert_hms_dec(dec):
    string = str(dec)
    if ':' in string:
        try:
            split = string.split(sep = ':')
            if float(split[0]) < 0:
                return float(split[0]) - float(split[1])/60 - float(split[2])/3600
            else:
                return float(split[0]) + float(split[1])/60 + float(split[2])/3600
        except:
            print('Cannot interpret ' + str(ra) + ' as a valid right ascension.')
            return -1
    else:
        return float(dec)
        
def convert_hms_ra(ra):
    string = str(ra)
    if ':' in string:
        try:
            split = string.split(sep = ':')
            return (float(split[0]) + float(split[1])/60 + float(split[2])/3600)*15.0
        except:
            print('Cannot interpret ' + str(ra) + ' as a valid right ascension.')
            return -1
    else:
        return float(ra)

def file_import(filename, delimiter = ','):
    if type(filename) == str and '.csv' in filename:
        data = pd.read_csv(filename, sep = delimiter)
        if len(data.columns) == 1:
            print('Cannot read objects from ' + filename + '. The default delimiter is set as a comma.')
            return -1
        objs = []
        for i in range(len(data['name'])):
            objs.append([data['ra'][i], data['dec'][i]])
        return objs
    
    elif type(filename) == list:
        return filename
    
    else:
        print('Cannot interpret objects.')
        return -1

def path(obj, obs):
    r, d = convert_hms_ra(obj[0]), convert_hms_dec(obj[1])
    a = obs.latitude
    
    angles = np.arange(0, 362, 1)  
    theta = angles-r
    
    return arcsin( sin(a)*sin(d) + cos(a)*cos(d)*cos(theta) )

def offset_path(obj, obs, mid_point):
    r, d = convert_hms_ra(obj[0]), convert_hms_dec(obj[1])
    a = obs.latitude
    
    angles = np.arange(mid_point-180, mid_point+180, 1)
    theta = angles-r
    
    return arcsin( sin(a)*sin(d) + cos(a)*cos(d)*cos(theta) )

def offset_path_pointing_restricted(obj, obs, mid_point, pointing_restrictions = []):
    angles = np.arange(mid_point-180, mid_point+180, 1)
    dirs, tau = calculate_pointing(obj, obs, angles)

    r, d = convert_hms_ra(obj[0]), convert_hms_dec(obj[1])
    a = obs.latitude
    theta = angles-r

    p = arcsin( sin(a)*sin(d) + cos(a)*cos(d)*cos(theta) )

    for i in range(len(dirs)):
        if dirs[i] in pointing_restrictions:
            p[i] = float('nan')

    return p

def solar_radec(date):
    for i, v in enumerate(vernal_zps):
        if date > v:
            vernal_index = i
            break
    for i, a in enumerate(autumnal_zps):
        if date > a:
            autumnal_index = i
            break
            
    if vernal_index > autumnal_index:
        autumnal_index += 1
        date_diff = date - vernal_zps[vernal_index]
        offset = 0
    else:
        vernal_index += 1
        date_diff = date - autumnal_zps[autumnal_index]
        offset = 180
    
    solar_ra = date_diff/abs(autumnal_zps[autumnal_index]-vernal_zps[vernal_index]) * 180 + offset
    solar_dec = (23.5*sin(solar_ra))
    
    return [solar_ra, solar_dec]

def lunar_radec(date):
    split = date.split(sep = '-')
    year, month, day = float(split[0]), float(split[1]), float(split[2])
    t = ts.utc(year, month, day)
    ra, dec, _ = earth.at(t).observe(moon).radec()
    
    return ra._degrees, dec.degrees

def twilights(date, obs):
    ra = np.arange(0, 362, 1)
    solar_path = path(solar_radec(date), obs)
    
    astronomical_twilight = [0, 0]
    nautical_twilight = [0, 0]
    civil_twilight = [0, 0]
    
    for i in range(len(ra) - 1):
        if solar_path[i] < -18.0 and solar_path[i+1] > -18.0:
            astronomical_twilight[1] = ra[i]
        if solar_path[i] > -18.0 and solar_path[i+1] < -18.0:
            astronomical_twilight[0] = ra[i]
            
    for i in range(len(ra) - 1):
        if solar_path[i] < -12.0 and solar_path[i+1] > -12.0:
            nautical_twilight[1] = ra[i]
        if solar_path[i] > -12.0 and solar_path[i+1] < -12.0:
            nautical_twilight[0] = ra[i]
            
    for i in range(len(ra) - 1):
        if solar_path[i] < -6.0 and solar_path[i+1] > -6.0:
            civil_twilight[1] = ra[i]
        if solar_path[i] > -6.0 and solar_path[i+1] < -6.0:
            civil_twilight[0] = ra[i]
            
    return astronomical_twilight, nautical_twilight, civil_twilight

def angle_to_moon(obj, lunar_coords):
    r_moon, d_moon = lunar_coords
    r_sn, d_sn = convert_hms_ra(obj[0]), convert_hms_dec(obj[1])
    angles = np.arange(0, 361, 1)
    theta_moon = angles - r_moon
    theta_sn = angles - r_sn

    return int(arccos( sin(d_moon)*sin(d_sn) + cos(d_moon)*cos(d_sn) * ( sin(theta_moon)*sin(theta_sn) + cos(theta_moon)*cos(theta_sn) ) )[0])

def calculate_pointing(obj, obs, angles):
    r, d = convert_hms_ra(obj[0]), convert_hms_dec(obj[1])
    a = obs.latitude
    
    theta = angles - r
    
    e = -cos(d)*sin(theta)
    n = (sin(d) * cos(a)) - (cos(d) * cos(theta) * sin(a))
    length = np.sqrt(e**2 + n**2)
    
    tau = arccos( n/length )
    
    dirs_e = ['N', 'NE', 'E', 'SE', 'S']
    dirs_w = ['N', 'NW', 'W', 'SW', 'S']
    
    ranges_e = {'N':[0.0, 22.5], 'NE':[22.5, 67.5], 'E':[67.5, 112.5], 'SE':[112.5, 157.5], 'S':[157.5, 180.0]}
    ranges_w = {'N':[0.0, 22.5], 'NW':[22.5, 67.5], 'W':[67.5, 112.5], 'SW':[112.5, 157.5], 'S':[157.5, 180.0]}
    
    dirs = []
    
    for i in range(len(e)):
        if e[i] >= 0.0:
            for index in dirs_e:
                if ranges_e[index][0]<=tau[i]<=ranges_e[index][1]:
                    dirs.append(index)
                    break
        else:
            for index in dirs_w:
                if ranges_w[index][0]<=tau[i]<=ranges_w[index][1]:
                    dirs.append(index)
                    break
                    
    return dirs, tau

def obs_plot(objs, obs_id, date = today_mjd, names = None, min_angle = None, max_angle = None, outfile = None, moon_angle = False, delimiter = ',', pointing_restrictions = []):
    obs = retrieve_observatory(obs_id)

    try:
        if type(date) == str:
            if '-' in date and 'T' not in date:
                date_ext = str(date) + 'T00:00.0'
                t = Time(date_ext, format = 'isot', scale = 'utc')
                date = float(t.mjd)
            else:
                t = Time(date, format = 'isot', scale = 'utc')
                date = float(t.mjd)


    except:
        print('The date parameter must be given in one of the following formats:\n\tmjd\n\tyyyy-mm-dd\n\tyyyy-mm-ddThh:mm:ss.ss')
        return -1

    date_ymd = Time(date, format = 'mjd', scale = 'utc').isot.split(sep = 'T')[0]

    lunar_coords = lunar_radec(date_ymd)

    if obs == -1:
        return

    if type(names) != list:
        if type(objs) == str and '.csv' in objs:
            try:
                names = list(pd.read_csv(objs, sep = delimiter)['name'])
            except:
                print('Cannot read object names from ' + objs + '. The default delimiter is set as a comma.')
    objs = file_import(objs, delimiter = delimiter)
    if objs == -1:
        return
    
    twi = twilights(date, obs)
    for i in range(len(twi)):   
        if twi[i][0] < twi[i][1]:
            twi[i][0] += 360
    
    midpoint_ras = np.arange(0, 362, 1)
    mid_point = np.mean([twi[0][0], twi[0][1]+360] )
    
    ra = np.arange(mid_point-180, mid_point+180, 1)
    
    fig = plt.figure(figsize = (15, 8), dpi = 200)
    ax = fig.add_subplot(111)
    
    for axis in ['top', 'bottom', 'right', 'left']:
        ax.spines[axis].set_linewidth(2)
    ax.xaxis.set_tick_params(width=2, length = 6)
    ax.yaxis.set_tick_params(width=2, length = 6)
    
    for j, obj in enumerate(objs):
        moon = angle_to_moon(obj, lunar_coords)
        if len(pointing_restrictions) == 0:
            obj_path = offset_path(obj, obs, mid_point)
        else:
            obj_path = offset_path_pointing_restricted(obj, obs, mid_point, pointing_restrictions = pointing_restrictions)

        if type(names) != list:
            ax.plot(ra, obj_path, lw = 4, label = str(obj), alpha = 0.6, linestyle = linestyles[int(j/10)])
        else:
            ax.plot(ra, obj_path, lw = 4, label = names[j], alpha = 0.6, linestyle = linestyles[int(j/10)])

        if moon_angle:
            for i in range(len(ra)):
                if (i+2*j)%40==0 and twi[2][0]-10<ra[i]<twi[2][1]+360+10 and 5<obj_path[i]<85:
                    ax.text(ra[i], obj_path[i], moon, va = 'center', ha = 'center')
    
    if min_angle != None:
        ax.axhline(min_angle, color = 'C3', linestyle = ':')
        ax.fill_between([mid_point-180, mid_point+180], 0, min_angle, color = 'C3', alpha = 0.1)
    if max_angle != None:
        ax.axhline(max_angle, color = 'C3', linestyle = ':')
        ax.fill_between([mid_point-180, mid_point+180], max_angle, 90, color = 'C3', alpha = 0.1)

    ax.plot(ra, offset_path(lunar_coords, obs, mid_point), lw = 3, color = 'k', linestyle = ':')

    ax.set_ylim(0, 90)
    ax.set_xlim(twi[2][0]-20, twi[2][1]+360+20)
    ax.legend(frameon = 0, ncol = 5, mode = 'expand', bbox_to_anchor=(0.0, -0.15, 1.0, 0.0))
    ax.set_xlabel('RA$_{peak}$')
    ax.set_ylabel(r'$\varphi$')
    ax.text(0.01, 1.01, obs.name + ' - ' + str(date_ymd), transform = ax.transAxes, va = 'bottom', ha = 'left', weight = 'bold')
    ax.text(0.99, 1.01, str(round(obs.latitude, 3)) + 'N, ' + str(round(obs.longitude, 3)) + 'E - Elevation: ' + str(obs.elevation) + 'm - ' + obs.country, transform = ax.transAxes, va = 'bottom', ha = 'right')
    if len(pointing_restrictions) != 0:
        line = 'Pointing restrictions: ' + pointing_restrictions[0]
        for pr in pointing_restrictions[1:]:
            line += ', ' + pr
        ax.text(0.99, 0.99, line, transform = ax.transAxes, va = 'top', ha = 'right')

    ax.axvline(twi[0][0], color = 'k', linestyle = ':')
    ax.fill_between([mid_point-180, twi[0][0]], 0, 90, color = 'k', alpha = 0.1)
    ax.axvline(twi[0][1]+360, color = 'k', linestyle = ':')
    ax.fill_between([twi[0][1]+360, mid_point+180], 0, 90, color = 'k', alpha = 0.1)
    
    ax.axvline(twi[1][0], color = 'k', linestyle = ':')
    ax.fill_between([mid_point-180, twi[1][0]], 0, 90, color = 'k', alpha = 0.1)
    ax.axvline(twi[1][1]+360, color = 'k', linestyle = ':')
    ax.fill_between([twi[1][1]+360, mid_point+180], 0, 90, color = 'k', alpha = 0.1)
    
    ax.axvline(twi[2][0], color = 'k', linestyle = ':')
    ax.fill_between([mid_point-180, twi[2][0]], 0, 90, color = 'k', alpha = 0.1)
    ax.axvline(twi[2][1]+360, color = 'k', linestyle = ':')
    ax.fill_between([twi[2][1]+360, mid_point+180], 0, 90, color = 'k', alpha = 0.1)
    
    
    plt.tight_layout()
    
    if outfile != None:
        plt.savefig(outfile, bbox_inches = 'tight')
        
    plt.show()

