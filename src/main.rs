use bio::io::fasta;
use cbor::Encoder;
use clap::Parser;
use crossbeam::{channel, thread};
use homopolymer_compress::{homopolymer_compress, homopolymer_compress_with_hodeco_map};
use log::{info, LevelFilter};
use simplelog::{ColorChoice, TermLogger, TerminalMode};
use std::fs::File;
use std::iter;
use std::path::PathBuf;

#[derive(Parser)]
struct Configuration {
    /// The input file.
    #[clap(index = 1, parse(from_os_str))]
    input: PathBuf,

    /// The output file. If not given, outputting to stdout.
    #[clap(index = 2, parse(from_os_str))]
    output: Option<PathBuf>,

    /// The file to output the map used to homopolymer decompress the output.
    #[clap(index = 3, parse(from_os_str))]
    hodeco_map_output: Option<PathBuf>,

    /// The number of compute threads to use for compressing.
    /// The program uses two extra threads for reading and writing the input and output files, which are not part of this number.
    /// It is likely that a very low number of threads is enough, since homopolymer compression is a very fast algorithm.
    #[clap(long, default_value = "1")]
    threads: usize,

    /// The size of the buffers between input and compute threads, and compute threads and output threads.
    #[clap(long, default_value = "32768")]
    buffer_size: usize,
}

fn initialise_logging() {
    TermLogger::init(
        LevelFilter::Debug,
        Default::default(),
        TerminalMode::Stderr,
        ColorChoice::Auto,
    )
    .unwrap();
    info!("Logging initialised successfully")
}

fn main() {
    let configuration = Configuration::parse();
    initialise_logging();

    let input = configuration.input;
    if let Some(extension) = input.extension() {
        if extension != "fasta" && extension != "fa" {
            panic!("Only fasta files supported at the moment, must end in .fa or .fasta, but ends in: {extension:?}");
        }
    } else {
        panic!("Only fasta files supported at the moment, must end in .fa or .fasta, but no extension found: {input:?}");
    }

    thread::scope(|scope| {
        let input_file =
            File::open(&input).unwrap_or_else(|error| panic!("Cannot open input file: {error:?}"));
        let (input_sender, input_receiver) = channel::bounded(configuration.buffer_size);
        scope
            .builder()
            .name("input_thread".to_string())
            .spawn(move |_| {
                for record in fasta::Reader::new(input_file).records() {
                    let record = record
                        .unwrap_or_else(|error| panic!("Cannot read fasta record: {error:?}"));
                    input_sender
                        .send(record)
                        .unwrap_or_else(|error| panic!("Cannot send fasta record: {error:?}"));
                }
            })
            .unwrap_or_else(|error| panic!("Cannot spawn input thread: {error:?}"));

        let (output_sender, output_receiver) =
            channel::bounded::<(String, Option<String>, (Vec<u8>, Option<Vec<usize>>))>(
                configuration.buffer_size,
            );
        let hodeco_map_output = configuration.hodeco_map_output.clone();
        if let Some(output) = configuration.output {
            let output_file = File::create(&output)
                .unwrap_or_else(|error| panic!("Cannot create output file: {error:?}"));
            scope
                .builder()
                .name("output_thread".to_string())
                .spawn(move |_| {
                    let mut writer = fasta::Writer::new(output_file);
                    let mut hodeco_mapping_writer = hodeco_map_output.as_ref().map(|path| {
                        Encoder::from_writer(File::create(path).unwrap_or_else(|error| {
                            panic!("Cannot create hodeco mapping output file: {error:?}")
                        }))
                    });
                    while let Ok((id, description, (sequence, hodeco_mapping))) =
                        output_receiver.recv()
                    {
                        writer
                            .write(&id, description.as_deref(), &sequence)
                            .unwrap_or_else(|error| panic!("Cannot write fasta record: {error:?}"));
                        if let Some(hodeco_mapping_writer) = hodeco_mapping_writer.as_mut() {
                            let hodeco_mapping = hodeco_mapping.unwrap_or_else(|| unreachable!());
                            hodeco_mapping_writer
                                .encode(iter::once((id, hodeco_mapping)))
                                .unwrap_or_else(|error| {
                                    panic!("Error writing hodeco mapping: {error:?}")
                                });
                        }
                    }
                })
                .unwrap_or_else(|error| panic!("Cannot spawn output thread: {error:?}"));
        } else {
            scope
                .builder()
                .name("output_thread".to_string())
                .spawn(move |_| {
                    let mut writer = fasta::Writer::new(std::io::stdout());
                    while let Ok((id, description, (sequence, hodeco_mapping))) =
                        output_receiver.recv()
                    {
                        assert!(
                            hodeco_mapping.is_none(),
                            "Found hodeco mapping even though no output file was specified."
                        );
                        writer
                            .write(&id, description.as_deref(), &sequence)
                            .unwrap_or_else(|error| panic!("Cannot write fasta record: {error:?}"));
                    }
                })
                .unwrap_or_else(|error| panic!("Cannot spawn output thread: {error:?}"));
        }

        for thread_id in 0..configuration.threads {
            let input_receiver = input_receiver.clone();
            let output_sender = output_sender.clone();
            let hodeco_map_output = configuration.hodeco_map_output.clone();
            scope
                .builder()
                .name(format!("compute_thread_{thread_id}"))
                .spawn(move |_| {
                    if hodeco_map_output.is_some() {
                        while let Ok(record) = input_receiver.recv() {
                            let (hoco_sequence, mut hodeco_mapping): (Vec<u8>, Vec<_>) =
                                homopolymer_compress_with_hodeco_map(record.seq().iter().cloned())
                                    .unzip();
                            hodeco_mapping.push(record.seq().len());
                            output_sender
                                .send((
                                    record.id().to_owned(),
                                    record.desc().map(str::to_owned),
                                    (hoco_sequence, Some(hodeco_mapping)),
                                ))
                                .unwrap_or_else(|error| {
                                    panic!("Cannot send fasta record: {error:?}")
                                });
                        }
                    } else {
                        while let Ok(record) = input_receiver.recv() {
                            let hoco_sequence: Vec<u8> =
                                homopolymer_compress(record.seq().iter().cloned()).collect();
                            output_sender
                                .send((
                                    record.id().to_owned(),
                                    record.desc().map(str::to_owned),
                                    (hoco_sequence, None),
                                ))
                                .unwrap_or_else(|error| {
                                    panic!("Cannot send fasta record: {error:?}")
                                });
                        }
                    }
                })
                .unwrap_or_else(|error| panic!("Cannot spawn compute thread: {error:?}"));
        }
    })
    .unwrap_or_else(|error| panic!("Error: {error:?}"));
}
