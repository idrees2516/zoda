# --- build stage -----------------------------------------------------------
FROM rust:1.83-slim AS builder
WORKDIR /build
COPY . .
RUN cargo build --release -p zoda-bench

# --- runtime stage ---------------------------------------------------------
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates && rm -rf /var/lib/apt/lists/*
COPY --from=builder /build/target/release/zoda-bench /usr/local/bin/zoda-bench
# Mount the mainnet trusted setup at /data/trusted_setup.txt
VOLUME ["/data"]
ENTRYPOINT ["zoda-bench"]
