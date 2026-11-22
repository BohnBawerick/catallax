// Nextest supplies the executable path after archive extraction.
macro_rules! bin_command {
    ($name:literal) => {
        std::process::Command::new(
            std::env::var_os(concat!("NEXTEST_BIN_EXE_", $name))
                .unwrap_or_else(|| env!(concat!("CARGO_BIN_EXE_", $name)).into()),
        )
    };
}

pub(crate) use bin_command;
