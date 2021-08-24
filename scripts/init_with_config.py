import os

def init_with_config(config: dict) -> dict:
    #########################################################
    # Set global path values and initialize exp object.
    #########################################################

    # The config.yaml must contain a name definition.
    if ('top_dir' not in config) and ('name' not in config):
        raise ValueError(
            f"Need name value in {configfile}. Prefer name & top_dir both set. top_dir=name when only name set.")

    # Set top_dir from name if not given.    
    if ('top_dir' not in config) and ('name' in config):
        config['top_dir'] = config['name']
        os.makedirs(config['top_dir'], exist_ok=True)

    # A star index is required.   
    if ('STAR_index' in config) and ('STAR index' not in config):
        config['STAR index'] = config['STAR_index']

    # Set some default path values as relative to the top_dir specified in config.yaml.
    _to = lambda x: config['top_dir'].rstrip('/') + f'/{x}'
    defaults = {
        'samples': _to('samples.txt'), 'beds': _to('beds'), 'fastq': _to('fastq'), 'sams': _to('sams'),
        'bigwig': _to('bigwig'), 'bedgraph': _to('bedgraph'), 'counts': _to('counts.txt'), 'STAR': 'STAR',
        'outs': _to('outs/'), 'counts': _to('outs/counts/'), 'data': _to('data'),
        'scheme': _to('samples.txt'),
    }

    for key in [_ for _ in defaults.keys() if (_ not in config)]:
        config[key] = defaults[key]

    # If a feature_gtf is not specified in the config.yaml, we give it a default
    # of assets/reference/gencode.v29.basic.annotation.gtf. The workflow will attempt
    # to download this file if it is not present. If the star index was made with a 
    # different annotation than gencode v29, some weird results might happen in
    # assigning reads to genes.
    if not os.path.exists(config['feature_gtf']):
        print("feature gtf not found: ", config['feature_gtf'])
        print("setting config['feature_gtf'] to match this default gtf: assets/reference/gencode.v29.basic.annotation.gtf")
        config['feature_gtf'] = "assets/reference/gencode.v29.basic.annotation.gtf"
        
    return config