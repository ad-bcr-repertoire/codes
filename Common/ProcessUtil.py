import configparser
import glob
import logging
import os
import platform
import subprocess

import pymysql

from .Enum import ChainType


def setup_database(sqlpath: str, dbname: str):
    """Create and setup a database schema to store g_target data."""

    logging.info('Started setup_database')

    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user='biteam',
                                 password='Antibody54321',
                                 port=3307,
                                 charset='utf8mb4',
                                 max_allowed_packet=1 * 1024 * 1024 * 1024,
                                 cursorclass=pymysql.cursors.DictCursor)

    # Create DB
    try:
        with connection.cursor() as cursor:
            cursor.execute("CREATE SCHEMA %s" % dbname)
            connection.commit()
    except Exception as e:
        logging.info('Pass database setup process.: %s', str(e))
        return
    finally:
        pass

    # Reconnect to the g_target database
    connection = pymysql.connect(host='localhost',
                                 db=dbname,
                                 user='biteam',
                                 password='Antibody54321',
                                 port=3307,
                                 charset='utf8mb4',
                                 max_allowed_packet=1 * 1024 * 1024 * 1024,
                                 cursorclass=pymysql.cursors.DictCursor)

    # Run setup sql scripts
    script_files = os.path.join(sqlpath, 'setup_*.sql')
    setup_sqls = glob.glob(script_files)
    for sql in setup_sqls:
        with open(sql, 'r') as f:
            sql_commands = f.read().split(';')
            for command in sql_commands:
                try:
                    with connection.cursor() as cursor:
                        cursor.execute(command)
                        connection.commit()
                except pymysql.err.InternalError:
                    pass
                finally:
                    pass

    logging.info('Finished setup_database')


def read_config(config_file: str) -> configparser.ConfigParser:
    directory = os.path.dirname(config_file)
    with open(config_file) as cf:
        config_str = cf.read()
        config_str = config_str.replace('$logfile', "r'%s'" % os.path.join(os.path.abspath(directory), 'log.log'))
    with open(config_file, 'w') as cf:
        cf.write(config_str)

    # read config
    config = configparser.ConfigParser()
    config.read(config_file)

    return config


def run_pear(forward: str, reverse: str, output: str, threads: int = 1) -> str:
    """Merge the paired-end sequence reads using PEAR - Paired-End reAd mergeR

    PEAR url: http://sco.h-its.org/exelixis/web/software/pear/
    """
    if platform.system() == 'Windows':
        pear_cmd = r'pear -f %s -r %s -o %s -j %d' % (forward, reverse, output, threads)
        process = subprocess.Popen(pear_cmd, stdout=subprocess.PIPE)
    elif platform.system() == 'Linux':
        script_directory = os.path.dirname(os.path.abspath(__file__))
        pear_cmd = ['pear', '-f', forward, '-r', reverse, '-o', output, '-j', str(threads)]
        process = subprocess.Popen(pear_cmd, stdout=subprocess.PIPE)
        process.communicate('Antibody54321\n')
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()

    return 1


def run_clustal_omega(input: str, output: str) -> str:
    if platform.system() == 'Windows':
        clustalo_cmd = r'clustalo -i %s -o %s --force' % (input, output)
        process = subprocess.Popen(clustalo_cmd, stdout=subprocess.DEVNULL)
    elif platform.system() == 'Linux':
        # script_directory = os.path.dirname(os.path.abspath(__file__))
        # clustalo_path = os.path.join(script_directory, 'clustalo-1.2.4-Ubuntu-x86_64')
        clustalo_path = r'/Tools/clustalo'
        # clustalo_cmd = ['sudo', '-S', clustalo_path, '-i', input, '-o', output, '--force']
        clustalo_cmd = [clustalo_path, '-i', input, '-o', output, '--force']
        process = subprocess.Popen(clustalo_cmd, stdout=subprocess.DEVNULL)
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()


import configparser
import glob
import logging
import os
import platform
import subprocess

import pymysql

from .Enum import ChainType


def setup_database(sqlpath: str, dbname: str):
    """Create and setup a database schema to store g_target data."""

    logging.info('Started setup_database')

    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user='biteam',
                                 password='Antibody54321',
                                 port=3307,
                                 charset='utf8mb4',
                                 max_allowed_packet=1 * 1024 * 1024 * 1024,
                                 cursorclass=pymysql.cursors.DictCursor)

    # Create DB
    try:
        with connection.cursor() as cursor:
            cursor.execute("CREATE SCHEMA %s" % dbname)
            connection.commit()
    except Exception as e:
        logging.info('Pass database setup process.: %s', str(e))
        return
    finally:
        pass

    # Reconnect to the g_target database
    connection = pymysql.connect(host='localhost',
                                 db=dbname,
                                 user='biteam',
                                 password='Antibody54321',
                                 port=3307,
                                 charset='utf8mb4',
                                 max_allowed_packet=1 * 1024 * 1024 * 1024,
                                 cursorclass=pymysql.cursors.DictCursor)

    # Run setup sql scripts
    script_files = os.path.join(sqlpath, 'setup_*.sql')
    setup_sqls = glob.glob(script_files)
    for sql in setup_sqls:
        with open(sql, 'r') as f:
            sql_commands = f.read().split(';')
            for command in sql_commands:
                try:
                    with connection.cursor() as cursor:
                        cursor.execute(command)
                        connection.commit()
                except pymysql.err.InternalError:
                    pass
                finally:
                    pass

    logging.info('Finished setup_database')


def read_config(config_file: str) -> configparser.ConfigParser:
    directory = os.path.dirname(config_file)
    with open(config_file) as cf:
        config_str = cf.read()
        config_str = config_str.replace('$logfile', "r'%s'" % os.path.join(os.path.abspath(directory), 'log.log'))
    with open(config_file, 'w') as cf:
        cf.write(config_str)

    # read config
    config = configparser.ConfigParser()
    config.read(config_file)

    return config


def run_pear(forward: str, reverse: str, output: str, threads: int = 1) -> str:
    """Merge the paired-end sequence reads using PEAR - Paired-End reAd mergeR

    PEAR url: http://sco.h-its.org/exelixis/web/software/pear/
    """
    if platform.system() == 'Windows':
        pear_cmd = r'pear -f %s -r %s -o %s -j %d' % (forward, reverse, output, threads)
        process = subprocess.Popen(pear_cmd, stdout=subprocess.PIPE)
    elif platform.system() == 'Linux':
        script_directory = os.path.dirname(os.path.abspath(__file__))
        pear_cmd = ['pear', '-f', forward, '-r', reverse, '-o', output, '-j', str(threads)]
        process = subprocess.Popen(pear_cmd, stdout=subprocess.PIPE)
        process.communicate('Antibody54321\n')
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()

    return 1


def run_clustal_omega(input: str, output: str) -> str:
    if platform.system() == 'Windows':
        clustalo_cmd = r'clustalo -i %s -o %s --force' % (input, output)
        process = subprocess.Popen(clustalo_cmd, stdout=subprocess.DEVNULL)
    elif platform.system() == 'Linux':
        script_directory = os.path.dirname(os.path.abspath(__file__))
        clustalo_path = os.path.join(script_directory, 'clustalo-1.2.4-Ubuntu-x86_64')
        clustalo_cmd = [clustalo_path, '-i', input, '-o', output, '--force']
        # process = subprocess.Popen(clustalo_cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        # process.stdin.write(b'aud2868ejr\n')
        # process.stdin.flush()
        # stdout, stderr = process.communicate(b'Antibody54321!')
        process = subprocess.Popen(clustalo_cmd, stdout=subprocess.DEVNULL)
        process.communicate()
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()


def run_igblast(chain_type: ChainType, query: str, out: str, domain_system: str = 'kabat', **kwargs):
    igblast_db_basepath = '/Tools/ncbi-igblast-1.8.0/bin/database'
    igblast_human_db_path = os.path.join(igblast_db_basepath, 'human')
    igblast_chicken_db_path = os.path.join(igblast_db_basepath, 'chicken')
    igblast_rabbit_db_path = os.path.join(igblast_db_basepath, 'rabbit')
    igblast_mouse_db_path = os.path.join(igblast_db_basepath, 'mouse')

    if chain_type in [ChainType.HUMAN_HEAVY, ChainType.HUMAN_LIGHT]:
        germline_db_V = os.path.join(igblast_human_db_path, r'imgt_human_v')
        germline_db_J = os.path.join(igblast_human_db_path, r'imgt_human_j')
        germline_db_D = os.path.join(igblast_human_db_path, r'imgt_human_d')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/human_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism human"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.HUMAN_BETA, ChainType.HUMAN_ALPHA, ChainType.HUMAN_DELTA, ChainType.HUMAN_GAMMA]:
        germline_db_V = os.path.join(igblast_human_db_path, r'imgt_human_trv_functional')
        germline_db_J = os.path.join(igblast_human_db_path, r'imgt_human_trj_functional')
        germline_db_D = os.path.join(igblast_human_db_path, r'imgt_human_trd_functional')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/human_gl.aux',
                                       query, out, domain_system) + \
              "-num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism human -ig_seqtype TCR -show_translation"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.CHICKEN_HEAVY, ChainType.CHICKEN_LIGHT]:
        germline_db_V = os.path.join(igblast_chicken_db_path, r'imgt_chicken_v')
        germline_db_J = os.path.join(igblast_chicken_db_path, r'imgt_chicken_j')
        germline_db_D = os.path.join(igblast_chicken_db_path, r'imgt_chicken_d')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/chicken_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism chicken"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.RABBIT_HEAVY, ChainType.RABBIT_KAPPA]:
        germline_db_V = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igv_whole')
        germline_db_J = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igj_whole')
        germline_db_D = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igd_whole')

        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/rabbit_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism rabbit"

        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.MOUSE_C57BL6_HEAVY]:
        germline_db_V = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_c57bl-6_igv_functional')
        germline_db_J = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_c57bl-6_igj_functional')
        germline_db_D = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_c57bl-6_igd_functional')

        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/mouse_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism mouse"

        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.MOUSE_BETA, ChainType.MOUSE_ALPHA]:
        germline_db_V = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trv_whole')
        germline_db_J = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trj_whole')
        germline_db_D = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trd_whole')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/mouse_gl.aux',
                                       query, out, domain_system) + \
              "-num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism mouse -ig_seqtype TCR -show_translation"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    if platform.system() == 'Windows':
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    elif platform.system() == 'Linux':
        igblast_path = os.path.join('/Tools', 'ncbi-igblast-1.8.0', 'bin', 'igblastn')
        # linux_cmd = ['sudo', '-S', igblast_path] + cmd.split(' ')[1:]
        linux_cmd = [igblast_path] + cmd.split(' ')[1:]
        process = subprocess.Popen(linux_cmd, stdout=subprocess.DEVNULL)
        # process.communicate('Antibody54321\n')
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()

def run_igblastp(chain_type: ChainType, query: str, out: str, domain_system: str = 'kabat', **kwargs):
    igblast_db_basepath = '/Tools/ncbi-igblast-1.8.0/bin/database'
    igblast_human_db_path = os.path.join(igblast_db_basepath, 'human')
    igblast_chicken_db_path = os.path.join(igblast_db_basepath, 'chicken')
    igblast_rabbit_db_path = os.path.join(igblast_db_basepath, 'rabbit')
    igblast_mouse_db_path = os.path.join(igblast_db_basepath, 'mouse')

    if chain_type in [ChainType.HUMAN_HEAVY, ChainType.HUMAN_LIGHT]:
        germline_db_V = os.path.join(igblast_human_db_path, r'imgt_human_v_prot')

        cmd = r"igblastp -germline_db_V %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, query, out, domain_system) + \
              "-num_alignments_V 1 -organism human"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.HUMAN_BETA, ChainType.HUMAN_ALPHA, ChainType.HUMAN_DELTA, ChainType.HUMAN_GAMMA]:
        germline_db_V = os.path.join(igblast_human_db_path, r'imgt_human_trv_functional')
        germline_db_J = os.path.join(igblast_human_db_path, r'imgt_human_trj_functional')
        germline_db_D = os.path.join(igblast_human_db_path, r'imgt_human_trd_functional')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/human_gl.aux',
                                       query, out, domain_system) + \
              "-num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism human -ig_seqtype TCR -show_translation"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.CHICKEN_HEAVY, ChainType.CHICKEN_LIGHT]:
        germline_db_V = os.path.join(igblast_chicken_db_path, r'imgt_chicken_v')
        germline_db_J = os.path.join(igblast_chicken_db_path, r'imgt_chicken_j')
        germline_db_D = os.path.join(igblast_chicken_db_path, r'imgt_chicken_d')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/chicken_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism chicken"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.RABBIT_HEAVY, ChainType.RABBIT_KAPPA]:
        germline_db_V = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igv_whole')
        germline_db_J = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igj_whole')
        germline_db_D = os.path.join(igblast_rabbit_db_path, r'imgt_oryctolagus_cuniculus_igd_whole')

        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/rabbit_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism rabbit"

        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.MOUSE_C57BL6_HEAVY]:
        germline_db_V = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_c57bl-6_igv_functional')
        germline_db_J = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_c57bl-6_igj_functional')
        germline_db_D = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_c57bl-6_igd_functional')

        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/mouse_gl.aux',
                                       query, out, domain_system) + \
              "-show_translation -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism mouse"

        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    elif chain_type in [ChainType.MOUSE_BETA, ChainType.MOUSE_ALPHA]:
        germline_db_V = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trv_whole')
        germline_db_J = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trj_whole')
        germline_db_D = os.path.join(igblast_mouse_db_path, r'imgt_mus_musculus_trd_whole')
        cmd = r"igblastn -germline_db_V %s -germline_db_J %s -germline_db_D %s -auxiliary_data %s -query %s -out %s " \
              r"-domain_system %s " % (germline_db_V, germline_db_J, germline_db_D, r'optional_file/mouse_gl.aux',
                                       query, out, domain_system) + \
              "-num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -organism mouse -ig_seqtype TCR -show_translation"
        for key, value in kwargs.items():
            cmd += ' -%s %s' % (key, str(value))

    if platform.system() == 'Windows':
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # elif platform.system() == 'Linux':
    #     igblast_path = os.path.join('/Tools', 'ncbi-igblast-1.8.0', 'bin', 'igblastp')
    #     linux_cmd = [igblast_path] + cmd.split(' ')[1:]
    #     print(linux_cmd)
    #     process = subprocess.Popen(linux_cmd, stdout=subprocess.DEVNULL)
    #     process.communicate('Antibody54321\n')
    elif platform.system() == 'Linux':
        igblast_path = os.path.join('/Tools', 'ncbi-igblast-1.8.0', 'bin', 'igblastn')
        # linux_cmd = ['sudo', '-S', igblast_path] + cmd.split(' ')[1:]
        linux_cmd = [igblast_path] + cmd.split(' ')[1:]
        process = subprocess.Popen(linux_cmd, stdout=subprocess.DEVNULL)
        # process.communicate('Antibody54321\n')
    else:
        logging.info('Not supporting OS')
        return None

    process.wait()

