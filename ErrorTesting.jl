# Error Testing Code

struct Animal
    name::String
    age::Int
end

struct Dog <: Animal
    breed::String
end

# Create an instance of the parent type
animal = Animal("Max", 3)

# Create an instance of the subtype
dog = Dog("Buddy", 5, "Golden Retriever")

# Access fields of the parent type
println(animal.name)  # Output: Max
println(animal.age)   # Output: 3

# Access fields of the subtype
println(dog.name)     # Output: Buddy
println(dog.age)      # Output: 5
println(dog.breed)    # Output: Golden Retriever
