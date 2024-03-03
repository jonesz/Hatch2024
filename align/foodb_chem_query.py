import mysql.connector

# Connect to the MySQL database
conn = mysql.connector.connect(
    host="localhost",
    user="username",
    password="password",
    database="your_database"
)
cursor = conn.cursor()

# Query to retrieve foods and their folic acid content
query = "SELECT food_name, folic_acid_amount FROM foods"

# Execute the query
cursor.execute(query)

# Fetch all results
foods = cursor.fetchall()

# Close the connection
conn.close()

# Print the table header
print("Food Name\t\tFolic Acid (mcg)")

# Print each food and its folic acid content
for food in foods:
    print(f"{food[0]}\t\t{food[1]} mcg")