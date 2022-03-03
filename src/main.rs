use bio::io::fasta;
use clap::Parser;
use crossbeam::{channel, thread};
use homopolymer_compress::homopolymer_compress;
use log::{info, LevelFilter};
use simplelog::{ColorChoice, TermLogger, TerminalMode};
use std::fs::File;
use std::path::PathBuf;

#[derive(Parser)]
struct Configuration {
    /// The input file.
    #[clap(index = 1, parse(from_os_str))]
    input: PathBuf,

    /// The output file. If not given, outputting to stdout.
    #[clap(index = 2, parse(from_os_str))]
    output: Option<PathBuf>,

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
            channel::bounded::<(String, Option<String>, Vec<u8>)>(configuration.buffer_size);
        if let Some(output) = configuration.output {
            let output_file = File::create(&output)
                .unwrap_or_else(|error| panic!("Cannot create output file: {error:?}"));
            scope
                .builder()
                .name("output_thread".to_string())
                .spawn(move |_| {
                    let mut writer = fasta::Writer::new(output_file);
                    while let Ok((id, description, sequence)) = output_receiver.recv() {
                        writer
                            .write(&id, description.as_deref(), &sequence)
                            .unwrap_or_else(|error| panic!("Cannot write fasta record: {error:?}"));
                    }
                })
                .unwrap_or_else(|error| panic!("Cannot spawn output thread: {error:?}"));
        } else {
            scope
                .builder()
                .name("output_thread".to_string())
                .spawn(move |_| {
                    let mut writer = fasta::Writer::new(std::io::stdout());
                    while let Ok((id, description, sequence)) = output_receiver.recv() {
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
            scope
                .builder()
                .name(format!("compute_thread_{thread_id}"))
                .spawn(move |_| {
                    while let Ok(record) = input_receiver.recv() {
                        let hoco_sequence: Vec<u8> =
                            homopolymer_compress(record.seq().iter().cloned()).collect();
                        output_sender
                            .send((
                                record.id().to_owned(),
                                record.desc().map(str::to_owned),
                                hoco_sequence,
                            ))
                            .unwrap_or_else(|error| panic!("Cannot send fasta record: {error:?}"));
                    }
                })
                .unwrap_or_else(|error| panic!("Cannot spawn compute thread: {error:?}"));
        }
    })
    .unwrap_or_else(|error| panic!("Error: {error:?}"));
}
