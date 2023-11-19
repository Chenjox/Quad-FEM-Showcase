

trait Chart {

  type Output;
  type Input;

  fn map_reference_to_actual(input_ref: Input) -> Output; 
}
