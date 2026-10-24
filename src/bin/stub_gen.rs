fn main() -> pyo3_stub_gen::Result<()> {
    fusion_neutron_utils::stub_info()?.generate()?;
    Ok(())
}
